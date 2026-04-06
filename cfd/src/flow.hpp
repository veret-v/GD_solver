#pragma once

#include <algorithm>
#include <cmath>

#include "mesh.hpp"

namespace cfd {

struct Primitive {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;
    double a = 0.0;
    double H = 0.0;
};

struct FlowParams {
    double Mach = 0.3;
    double alpha = 0.0;
    double gamma = 1.4;
    double rho_inf = 1.0;
    double p_inf = 1.0 / 1.4;
    double a_inf = 0.0;
    double u_inf = 0.0;
    double v_inf = 0.0;

    void init()
    {
        a_inf = std::sqrt(gamma * p_inf / rho_inf);
        const double V = Mach * a_inf;
        u_inf = V * std::cos(alpha);
        v_inf = V * std::sin(alpha);
    }

    Vec4 Uinf() const;
};

inline FlowParams FLOW{};

inline double pressure(const Vec4 &U, double g = FLOW.gamma)
{
    const double rho = std::max(U[0], 1e-12);
    const double kinetic = 0.5 * (U[1] * U[1] + U[2] * U[2]) / rho;
    return (g - 1.0) * (U[3] - kinetic);
}

inline double sound_speed(const Vec4 &U, double g = FLOW.gamma)
{
    return std::sqrt(g * std::max(pressure(U, g), 1e-12) / std::max(U[0], 1e-12));
}

inline Vec4 primitive_to_conservative(double rho, double u, double v, double p, double g = FLOW.gamma)
{
    const double E = p / (g - 1.0) + 0.5 * rho * (u * u + v * v);
    return make_vec4(rho, rho * u, rho * v, E);
}

inline Primitive conservative_to_primitive(const Vec4 &U, double g = FLOW.gamma)
{
    Primitive W;
    W.rho = std::max(U[0], 1e-12);
    W.u = U[1] / W.rho;
    W.v = U[2] / W.rho;
    W.p = std::max(pressure(U, g), 1e-12);
    W.a = std::sqrt(g * W.p / W.rho);
    W.H = (U[3] + W.p) / W.rho;
    return W;
}

inline void enforce_physical(Vec4 &U, double g = FLOW.gamma)
{
    U[0] = std::max(U[0], 1e-10);
    const double u = U[1] / U[0];
    const double v = U[2] / U[0];
    const double kinetic = 0.5 * U[0] * (u * u + v * v);
    const double Emin = kinetic + 1e-10 / (g - 1.0);
    U[3] = std::max(U[3], Emin);
}

inline Vec4 FlowParams::Uinf() const
{
    return primitive_to_conservative(rho_inf, u_inf, v_inf, p_inf, gamma);
}

} // namespace cfd
