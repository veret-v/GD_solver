#pragma once

#include <fstream>
#include <stdexcept>
#include <string>

#include "flow.hpp"

namespace cfd {

inline void init_freestream(Mesh &m)
{
    const Vec4 Ui = FLOW.Uinf();
    for(auto &c : m.cells) c.U = Ui;
}

inline void init_potential(Mesh &m, double cx, double cy, double R, double g = FLOW.gamma)
{
    const double U0 = std::sqrt(FLOW.u_inf * FLOW.u_inf + FLOW.v_inf * FLOW.v_inf);
    for(auto &c : m.cells) {
        const double dx = c.cx - cx;
        const double dy = c.cy - cy;
        const double r2 = std::max(dx * dx + dy * dy, 1.05 * R * R);
        const double R2r2 = R * R / r2;
        const double u = FLOW.u_inf * (1.0 - R2r2 * (dx * dx - dy * dy) / r2) -
                         2.0 * FLOW.v_inf * R2r2 * dx * dy / r2;
        const double v = FLOW.v_inf * (1.0 + R2r2 * (dx * dx - dy * dy) / r2) -
                         2.0 * FLOW.u_inf * R2r2 * dx * dy / r2;
        const double p = std::max(FLOW.p_inf +
            0.5 * FLOW.rho_inf * (U0 * U0 - u * u - v * v), 1e-6);
        c.U = primitive_to_conservative(FLOW.rho_inf, u, v, p, g);
    }
}

inline void load_restart(Mesh &m, const std::string &filename)
{
    std::ifstream f(filename, std::ios::binary);
    if(!f) throw std::runtime_error("Failed to open restart file: " + filename);

    int nc = 0;
    f.read(reinterpret_cast<char *>(&nc), sizeof(int));
    if(nc != m.nc()) throw std::runtime_error("Restart cell count mismatch");

    for(auto &c : m.cells) {
        f.read(reinterpret_cast<char *>(c.U.data()), static_cast<std::streamsize>(4 * sizeof(double)));
        if(!f) throw std::runtime_error("Restart file is truncated");
        enforce_physical(c.U);
    }
}

} // namespace cfd
