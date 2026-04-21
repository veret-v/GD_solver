#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>

#include "flow.hpp"

namespace cfd {

constexpr double kRadToDeg = 57.2957795130823208768;

inline void write_vtk(const Mesh &m, const std::string &filename, double g = FLOW.gamma)
{
    std::ofstream f(filename);
    if(!f) throw std::runtime_error("Failed to open VTK file: " + filename);

    f << std::scientific << std::setprecision(17);
    f << "# vtk DataFile Version 3.0\n";
    f << "CFD\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";

    f << "POINTS " << m.nodes.size() << " double\n";
    for(const Node &n : m.nodes) f << n.x << " " << n.y << " 0\n";

    const int nc = m.n_owned > 0 ? m.n_owned : m.nc();
    f << "CELLS " << nc << " " << nc * 4 << "\n";
    for(int ci = 0; ci < nc; ++ci) {
        const Cell &c = m.cells[ci];
        f << "3";
        for(int id : c.node_ids) f << " " << id;
        f << "\n";
    }

    f << "CELL_TYPES " << nc << "\n";
    for(int i = 0; i < nc; ++i) f << "5\n";

    f << "CELL_DATA " << nc << "\n";

    f << "SCALARS density double 1\nLOOKUP_TABLE default\n";
    for(int ci = 0; ci < nc; ++ci) f << m.cells[ci].U[0] << "\n";

    f << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for(int ci = 0; ci < nc; ++ci) f << conservative_to_primitive(m.cells[ci].U, g).p << "\n";

    f << "SCALARS mach double 1\nLOOKUP_TABLE default\n";
    for(int ci = 0; ci < nc; ++ci) {
        const Cell &c = m.cells[ci];
        const Primitive W = conservative_to_primitive(c.U, g);
        f << std::sqrt(W.u * W.u + W.v * W.v) / W.a << "\n";
    }

    f << "VECTORS velocity double\n";
    for(int ci = 0; ci < nc; ++ci) {
        const Cell &c = m.cells[ci];
        const Primitive W = conservative_to_primitive(c.U, g);
        f << W.u << " " << W.v << " 0\n";
    }
}


inline std::string snap_name(const std::string &results_dir, int snap)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ostringstream oss;
    oss << results_dir << "/out_" << std::setfill('0') << std::setw(4) << snap 
        << "_r" << rank << ".vtk"; // Добавили _r0, _r1...
    return oss.str();
}

inline void write_snap(const Mesh &m, const std::string &results_dir, int snap, double g = FLOW.gamma)
{
    
    std::filesystem::create_directories(results_dir);
    write_vtk(m, snap_name(results_dir, snap), g);
}

inline void save_restart(const Mesh &m, const std::string &filename, int rank)
{
    // Формируем имя файла с суффиксом ранга
    
    std::string ranked_name = filename + "_r" + std::to_string(rank) + ".bin";
    std::ofstream f(ranked_name, std::ios::binary);
    
    if(!f) throw std::runtime_error("Failed to open restart file: " + ranked_name);

    const int nc = m.nc();
    f.write(reinterpret_cast<const char *>(&nc), sizeof(int));
    for(const Cell &c : m.cells) {
        f.write(reinterpret_cast<const char *>(c.U.data()),
                static_cast<std::streamsize>(4 * sizeof(double)));
    }
}

inline void write_cp(const Mesh &m,
                     const std::string &filename,
                     int rank, // добавили rank
                     double center_x,
                     double center_y,
                     double rho_inf,
                     double V_inf,
                     double p_inf,
                     double g = FLOW.gamma)
{
    
    std::string ranked_name = filename + "_r" + std::to_string(rank) + ".dat";
    std::ofstream f(ranked_name); // создаем файл для конкретного ранга
    const double q = 0.5 * rho_inf * V_inf * V_inf;
    std::vector<std::pair<double, double>> data;
    data.reserve(m.faces.size());

    for(const Face &f : m.faces) {
        if(f.bc != Face::BC::Wall) continue;
        const Primitive W = conservative_to_primitive(m.cells[f.left].U, g);
        const double theta = std::atan2(f.my - center_y, f.mx - center_x) * kRadToDeg;
        const double cp = (W.p - p_inf) / q;
        data.emplace_back(theta, cp);
    }

    std::sort(data.begin(), data.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    
    if(!f) throw std::runtime_error("Failed to open Cp file: " + filename);

    f << "# theta_deg Cp\n";
    for(const auto &[theta, cp] : data) f << theta << " " << cp << "\n";
}

} // namespace cfd
