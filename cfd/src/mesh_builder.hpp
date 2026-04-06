#pragma once

#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cmath>

#include <gmsh.h>

#include "mesh.hpp"

namespace cfd {

inline std::string to_lower(std::string s)
{
    for(char &ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

inline Face::BC bc_from_name(const std::string &name)
{
    const std::string low = to_lower(name);
    if(low == "wall" || low == "airfoil" || low == "cylinder") return Face::BC::Wall;
    if(low == "inflow" || low == "inlet") return Face::BC::Inflow;
    if(low == "outflow" || low == "outlet") return Face::BC::Outflow;
    if(low == "farfield" || low == "far_field") return Face::BC::Farfield;
    if(low == "symmetry") return Face::BC::Symmetry;
    return Face::BC::Interior;
}

inline void cell_geometry(Mesh &m)
{
    for(auto &c : m.cells) {
        if(c.node_ids.size() < 3) throw std::runtime_error("Cell has fewer than 3 nodes");
        double A = 0.0;
        double cx = 0.0;
        double cy = 0.0;
        const int n = static_cast<int>(c.node_ids.size());
        for(int k = 0; k < n; ++k) {
            const Node &A0 = m.nodes[c.node_ids[k]];
            const Node &B0 = m.nodes[c.node_ids[(k + 1) % n]];
            const double cr = A0.x * B0.y - B0.x * A0.y;
            A += cr;
            cx += (A0.x + B0.x) * cr;
            cy += (A0.y + B0.y) * cr;
        }
        A *= 0.5;
        if(std::abs(A) < 1e-14) throw std::runtime_error("Degenerate cell area");
        c.vol = std::abs(A);
        c.cx = cx / (6.0 * A);
        c.cy = cy / (6.0 * A);
    }
}

inline void build_faces(Mesh &m)
{
    std::map<std::pair<int, int>, std::vector<int>> edge_map;
    for(int ci = 0; ci < m.nc(); ++ci) {
        auto &cell = m.cells[ci];
        cell.face_ids.clear();
        const int n = static_cast<int>(cell.node_ids.size());
        for(int k = 0; k < n; ++k) {
            const int a = cell.node_ids[k];
            const int b = cell.node_ids[(k + 1) % n];
            edge_map[{std::min(a, b), std::max(a, b)}].push_back(ci);
        }
    }

    m.faces.clear();
    m.faces.reserve(edge_map.size());

    for(const auto &[key, owners] : edge_map) {
        if(owners.empty() || owners.size() > 2) {
            throw std::runtime_error("Unsupported non-manifold edge in build_faces()");
        }

        Face f;
        f.node0 = key.first;
        f.node1 = key.second;
        const Node &A = m.nodes[f.node0];
        const Node &B = m.nodes[f.node1];

        const double dx = B.x - A.x;
        const double dy = B.y - A.y;
        f.length = std::sqrt(dx * dx + dy * dy);
        f.mx = 0.5 * (A.x + B.x);
        f.my = 0.5 * (A.y + B.y);
        f.nx = dy / f.length;
        f.ny = -dx / f.length;

        if(owners.size() == 2) {
            f.left = owners[0];
            f.right = owners[1];
            const double dcx = m.cells[f.right].cx - m.cells[f.left].cx;
            const double dcy = m.cells[f.right].cy - m.cells[f.left].cy;
            if(f.nx * dcx + f.ny * dcy < 0.0) {
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
            f.bc = Face::BC::Interior;
        } else {
            f.left = owners[0];
            f.right = -1;
            const double dot = f.nx * (f.mx - m.cells[f.left].cx) +
                               f.ny * (f.my - m.cells[f.left].cy);
            if(dot < 0.0) {
                f.nx = -f.nx;
                f.ny = -f.ny;
            }
            f.bc = Face::BC::Wall;
        }

        const int fid = static_cast<int>(m.faces.size());
        m.cells[f.left].face_ids.push_back(fid);
        if(f.right >= 0) m.cells[f.right].face_ids.push_back(fid);
        m.faces.push_back(f);
    }
}

class GmshSession {
public:
    GmshSession()
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", 1);
    }

    ~GmshSession()
    {
        try {
            gmsh::clear();
            gmsh::finalize();
        } catch(...) {
        }
    }

    GmshSession(const GmshSession &) = delete;
    GmshSession &operator=(const GmshSession &) = delete;
};

inline Mesh load_gmsh(const std::string &filename)
{
    GmshSession session;
    gmsh::open(filename);

    Mesh m;

    std::vector<std::size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> parametric;
    gmsh::model::mesh::getNodes(node_tags, coords, parametric);

    m.nodes.resize(node_tags.size());
    std::unordered_map<std::size_t, int> node_to_local;
    node_to_local.reserve(node_tags.size());

    for(std::size_t i = 0; i < node_tags.size(); ++i) {
        m.nodes[i] = Node{coords[3 * i + 0], coords[3 * i + 1]};
        node_to_local[node_tags[i]] = static_cast<int>(i);
    }

    std::map<std::pair<int, int>, Face::BC> edge_bc;
    std::vector<std::pair<int, int>> phys_groups;
    gmsh::model::getPhysicalGroups(phys_groups, 1);
    for(const auto &[dim, phys_tag] : phys_groups) {
        std::string name;
        gmsh::model::getPhysicalName(dim, phys_tag, name);
        const Face::BC bc = bc_from_name(name);

        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(dim, phys_tag, entities);
        for(const int entity_tag : entities) {
            std::vector<int> elem_types;
            std::vector<std::vector<std::size_t>> elem_tags;
            std::vector<std::vector<std::size_t>> elem_nodes;
            gmsh::model::mesh::getElements(elem_types, elem_tags, elem_nodes, 1, entity_tag);

            for(std::size_t it = 0; it < elem_types.size(); ++it) {
                const int type = elem_types[it];
                const std::size_t nodes_per_elem = (type == 1) ? 2u : (type == 8 ? 3u : 0u);
                if(nodes_per_elem == 0) continue;
                const auto &nodes = elem_nodes[it];
                const std::size_t ne = nodes.size() / nodes_per_elem;
                for(std::size_t e = 0; e < ne; ++e) {
                    const std::size_t gmsh_n0 = nodes[e * nodes_per_elem + 0];
                    const std::size_t gmsh_n1 = nodes[e * nodes_per_elem + (nodes_per_elem - 1)];
                    const int n0 = node_to_local.at(gmsh_n0);
                    const int n1 = node_to_local.at(gmsh_n1);
                    edge_bc[{std::min(n0, n1), std::max(n0, n1)}] = bc;
                }
            }
        }
    }

    std::vector<int> elem_types;
    std::vector<std::vector<std::size_t>> elem_tags;
    std::vector<std::vector<std::size_t>> elem_nodes;
    gmsh::model::mesh::getElements(elem_types, elem_tags, elem_nodes, 2, -1);
    for(std::size_t it = 0; it < elem_types.size(); ++it) {
        const int type = elem_types[it];
        if(type != 2) continue;
        const auto &nodes = elem_nodes[it];
        const std::size_t ne = nodes.size() / 3u;
        for(std::size_t e = 0; e < ne; ++e) {
            Cell c;
            c.node_ids = {
                node_to_local.at(nodes[3 * e + 0]),
                node_to_local.at(nodes[3 * e + 1]),
                node_to_local.at(nodes[3 * e + 2])
            };
            m.cells.push_back(c);
        }
    }

    if(m.cells.empty()) throw std::runtime_error("No triangular cells found in Gmsh file");

    cell_geometry(m);
    build_faces(m);

    for(auto &f : m.faces) {
        if(!f.is_boundary()) continue;
        const auto key = std::make_pair(std::min(f.node0, f.node1), std::max(f.node0, f.node1));
        const auto it = edge_bc.find(key);
        if(it != edge_bc.end()) f.bc = it->second;
    }

    return m;
}

} // namespace cfd
