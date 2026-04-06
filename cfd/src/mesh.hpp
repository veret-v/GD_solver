#pragma once

#include <array>
#include <vector>

namespace cfd {

using Vec4 = std::array<double, 4>;

inline Vec4 make_vec4(double a = 0.0, double b = 0.0, double c = 0.0, double d = 0.0)
{
    return Vec4{a, b, c, d};
}

inline Vec4 &operator+=(Vec4 &lhs, const Vec4 &rhs)
{
    for(int i = 0; i < 4; ++i) lhs[i] += rhs[i];
    return lhs;
}

inline Vec4 &operator-=(Vec4 &lhs, const Vec4 &rhs)
{
    for(int i = 0; i < 4; ++i) lhs[i] -= rhs[i];
    return lhs;
}

inline Vec4 &operator*=(Vec4 &lhs, double s)
{
    for(double &v : lhs) v *= s;
    return lhs;
}

inline Vec4 &operator/=(Vec4 &lhs, double s)
{
    for(double &v : lhs) v /= s;
    return lhs;
}

inline Vec4 operator+(Vec4 lhs, const Vec4 &rhs) { return lhs += rhs; }
inline Vec4 operator-(Vec4 lhs, const Vec4 &rhs) { return lhs -= rhs; }
inline Vec4 operator*(Vec4 lhs, double s) { return lhs *= s; }
inline Vec4 operator*(double s, Vec4 rhs) { return rhs *= s; }
inline Vec4 operator/(Vec4 lhs, double s) { return lhs /= s; }

struct Node {
    double x = 0.0;
    double y = 0.0;
};

struct Face {
    enum class BC {
        Interior,
        Wall,
        Inflow,
        Outflow,
        Farfield,
        Symmetry
    };

    int left = -1;
    int right = -1;
    int node0 = -1;
    int node1 = -1;
    double nx = 0.0;
    double ny = 0.0;
    double length = 0.0;
    double mx = 0.0;
    double my = 0.0;
    BC bc = BC::Interior;

    bool is_boundary() const { return right < 0; }
};

struct Cell {
    std::vector<int> face_ids;
    std::vector<int> node_ids;
    double vol = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    Vec4 U = make_vec4();
    Vec4 res = make_vec4();
};

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

    int nc() const { return static_cast<int>(cells.size()); }
    int nf() const { return static_cast<int>(faces.size()); }

    void zero_res()
    {
        for(auto &c : cells) c.res = make_vec4();
    }
};

inline int other_cell(const Face &f, int cell_id)
{
    if(f.left == cell_id) return f.right;
    if(f.right == cell_id) return f.left;
    return -1;
}

} // namespace cfd
