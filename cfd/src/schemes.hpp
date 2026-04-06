#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "bc.hpp"
#include "riemann.hpp"

namespace cfd {

enum class SpatialScheme {
    Godunov
};

enum class TimeScheme {
    Euler,
    RK2
};

struct Aero {
    double Cd = 0.0;
    double Cl = 0.0;
};

inline void residuals_godunov(Mesh &m, double g = FLOW.gamma)
{
    m.zero_res();
    for(const Face &f : m.faces) {
        const Vec4 &UL = m.cells[f.left].U;
        const Vec4 UR = f.is_boundary()
            ? ghost_state(UL, f.nx, f.ny, f.bc, g)
            : m.cells[f.right].U;
        const Vec4 flux = hllc(UL, UR, f.nx, f.ny, g);

        m.cells[f.left].res -= flux * f.length;
        if(f.right >= 0) m.cells[f.right].res += flux * f.length;
    }
}

inline void compute_residuals(Mesh &m, SpatialScheme scheme, double g = FLOW.gamma)
{
    switch(scheme) {
    case SpatialScheme::Godunov:
        residuals_godunov(m, g);
        break;
    }
}

inline void update_cells(Mesh &m, double dt, double g = FLOW.gamma)
{
    for(Cell &c : m.cells) {
        c.U += (dt / c.vol) * c.res;
        enforce_physical(c.U, g);
    }
}

inline void step_euler(Mesh &m, double dt, SpatialScheme scheme, double g = FLOW.gamma)
{
    compute_residuals(m, scheme, g);
    update_cells(m, dt, g);
}

inline void step_rk2(Mesh &m, double dt, SpatialScheme scheme, double g = FLOW.gamma)
{
    Mesh stage = m;

    compute_residuals(m, scheme, g);
    for(int i = 0; i < m.nc(); ++i) {
        stage.cells[i].U = m.cells[i].U + (dt / m.cells[i].vol) * m.cells[i].res;
        enforce_physical(stage.cells[i].U, g);
    }

    compute_residuals(stage, scheme, g);
    for(int i = 0; i < m.nc(); ++i) {
        const Vec4 U1 = stage.cells[i].U + (dt / stage.cells[i].vol) * stage.cells[i].res;
        m.cells[i].U = 0.5 * (m.cells[i].U + U1);
        enforce_physical(m.cells[i].U, g);
    }
}

inline void step(Mesh &m, double dt, SpatialScheme scheme, TimeScheme time_scheme, double g = FLOW.gamma)
{
    switch(time_scheme) {
    case TimeScheme::Euler:
        step_euler(m, dt, scheme, g);
        break;
    case TimeScheme::RK2:
        step_rk2(m, dt, scheme, g);
        break;
    }
}

inline double compute_dt(const Mesh &m, double CFL, double g = FLOW.gamma)
{
    double dt = 1e100;
    for(int ci = 0; ci < m.nc(); ++ci) {
        const Cell &c = m.cells[ci];
        const Primitive W = conservative_to_primitive(c.U, g);

        double sigma = 0.0;
        for(const int fid : c.face_ids) {
            const Face &f = m.faces[fid];
            const double lam = std::abs(W.u * f.nx + W.v * f.ny) + W.a;
            sigma += lam * f.length;
        }

        dt = std::min(dt, CFL * c.vol / (sigma + 1e-12));
    }
    if(!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::runtime_error("compute_dt() produced invalid dt");
    }
    return dt;
}

inline Aero compute_aero(const Mesh &m, double rho_inf, double V_inf, double D_ref, double g = FLOW.gamma)
{
    double Fx = 0.0;
    double Fy = 0.0;
    const double q = 0.5 * rho_inf * V_inf * V_inf;

    for(const Face &f : m.faces) {
        if(f.bc != Face::BC::Wall) continue;
        const Primitive W = conservative_to_primitive(m.cells[f.left].U, g);
        Fx -= W.p * f.nx * f.length;
        Fy -= W.p * f.ny * f.length;
    }

    return Aero{Fx / (q * D_ref), Fy / (q * D_ref)};
}

} // namespace cfd
