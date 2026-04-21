#pragma once

#include <algorithm>
#include <cmath>

#include "flow.hpp"

namespace cfd {

inline Vec4 bc_wall(const Vec4 &UL, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive W = conservative_to_primitive(UL, g);
    const double un = W.u * nx + W.v * ny;
    const double ug = W.u - 2.0 * un * nx;
    const double vg = W.v - 2.0 * un * ny;
    return primitive_to_conservative(W.rho, ug, vg, W.p, g);
}

inline Vec4 bc_inflow(double g = FLOW.gamma)
{
    (void)g;
    return FLOW.Uinf();
}

inline Vec4 bc_outflow(const Vec4 &UL)
{
    return UL;
}

inline Vec4 bc_farfield(const Vec4 &UL, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive L = conservative_to_primitive(UL, g);
    const double tx = -ny;
    const double ty = nx;

    const double unL = L.u * nx + L.v * ny;
    const double utL = L.u * tx + L.v * ty;

    const double unInf = FLOW.u_inf * nx + FLOW.v_inf * ny;
    const double utInf = FLOW.u_inf * tx + FLOW.v_inf * ty;

    const double Rp = unL + 2.0 * L.a / (g - 1.0);
    const double Rm = unInf - 2.0 * FLOW.a_inf / (g - 1.0);

    const double unG = 0.5 * (Rp + Rm);
    const double aG = std::max(0.25 * (g - 1.0) * (Rp - Rm), 1e-8);
    const bool inflow_dominated = (unL < 0.0);
    const double sG = inflow_dominated
        ? FLOW.p_inf / std::pow(FLOW.rho_inf, g)
        : L.p / std::pow(L.rho, g);

    const double rhoG = std::pow((aG * aG) / (g * sG), 1.0 / (g - 1.0));
    const double pG = sG * std::pow(rhoG, g);
    const double utG = inflow_dominated ? utInf : utL;

    const double uG = unG * nx + utG * tx;
    const double vG = unG * ny + utG * ty;
    return primitive_to_conservative(rhoG, uG, vG, pG, g);
}

// CHECK: GHOST_STATE
inline Vec4 ghost_state(const Vec4 &UL, double nx, double ny, Face::BC bc, double g = FLOW.gamma)
{
    switch(bc) {
    case Face::BC::MPIBound: return UL;
    case Face::BC::Interior:
        return UL;
    case Face::BC::Wall:
    case Face::BC::Symmetry:
        return bc_wall(UL, nx, ny, g);
    case Face::BC::Inflow:
        return bc_inflow(g);
    case Face::BC::Outflow:
        return bc_outflow(UL);
    case Face::BC::Farfield:
        return bc_farfield(UL, nx, ny, g);
    }
    
    return UL;
}

} // namespace cfd
