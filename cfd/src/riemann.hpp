#pragma once

#include <algorithm>
#include <cmath>

#include "flow.hpp"

namespace cfd {

inline Vec4 flux_normal(const Vec4 &U, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive W = conservative_to_primitive(U, g);
    const double un = W.u * nx + W.v * ny;
    return make_vec4(
        W.rho * un,
        W.rho * W.u * un + W.p * nx,
        W.rho * W.v * un + W.p * ny,
        (U[3] + W.p) * un
    );
}

inline Vec4 hlle(const Vec4 &UL, const Vec4 &UR, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive L = conservative_to_primitive(UL, g);
    const Primitive R = conservative_to_primitive(UR, g);

    const double unL = L.u * nx + L.v * ny;
    const double unR = R.u * nx + R.v * ny;
    const double SL = std::min(unL - L.a, unR - R.a);
    const double SR = std::max(unL + L.a, unR + R.a);

    const Vec4 FL = flux_normal(UL, nx, ny, g);
    const Vec4 FR = flux_normal(UR, nx, ny, g);

    if(0.0 <= SL) return FL;
    if(SR <= 0.0) return FR;

    return (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL + 1e-300);
}

inline Vec4 star_state(const Vec4 &U, double S, double Sstar, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive W = conservative_to_primitive(U, g);
    const double un = W.u * nx + W.v * ny;
    const double coeff = W.rho * (S - un) / (S - Sstar);
    const double us = W.u + (Sstar - un) * nx;
    const double vs = W.v + (Sstar - un) * ny;
    const double Es = U[3] / W.rho +
                      (Sstar - un) * (Sstar + W.p / (W.rho * (S - un)));
    return make_vec4(coeff, coeff * us, coeff * vs, coeff * Es);
}

inline Vec4 hllc(const Vec4 &UL, const Vec4 &UR, double nx, double ny, double g = FLOW.gamma)
{
    const Primitive L = conservative_to_primitive(UL, g);
    const Primitive R = conservative_to_primitive(UR, g);

    const double unL = L.u * nx + L.v * ny;
    const double unR = R.u * nx + R.v * ny;

    const double sqL = std::sqrt(L.rho);
    const double sqR = std::sqrt(R.rho);
    const double inv = 1.0 / (sqL + sqR);

    const double uRoe = (sqL * L.u + sqR * R.u) * inv;
    const double vRoe = (sqL * L.v + sqR * R.v) * inv;
    const double HRoe = (sqL * L.H + sqR * R.H) * inv;
    const double unRoe = uRoe * nx + vRoe * ny;
    const double aRoe = std::sqrt(std::max((g - 1.0) * (HRoe - 0.5 * (uRoe * uRoe + vRoe * vRoe)), 1e-12));

    const double SL = std::min({unL - L.a, unRoe - aRoe});
    const double SR = std::max({unR + R.a, unRoe + aRoe});

    const double denom = L.rho * (SL - unL) - R.rho * (SR - unR);
    //if(std::abs(denom) < 1e-12 || !(SL < SR)) {
    //    return hlle(UL, UR, nx, ny, g);
    //}

    const double Sstar = (R.p - L.p +
        L.rho * unL * (SL - unL) -
        R.rho * unR * (SR - unR)) / (denom + 1e-300);

    const Vec4 FL = flux_normal(UL, nx, ny, g);
    const Vec4 FR = flux_normal(UR, nx, ny, g);

    if(0.0 <= SL) return FL;
    if(SR <= 0.0) return FR;

    if(0.0 <= Sstar) {
        const Vec4 UstarL = star_state(UL, SL, Sstar, nx, ny, g);
        //if(!(UstarL[0] > 0.0) || !std::isfinite(UstarL[3])) {
         //   return hlle(UL, UR, nx, ny, g);
        //}
        return FL + SL * (UstarL - UL);
    }

    const Vec4 UstarR = star_state(UR, SR, Sstar, nx, ny, g);
    if(!(UstarR[0] > 0.0) || !std::isfinite(UstarR[3])) {
        return hlle(UL, UR, nx, ny, g);
    }
    return FR + SR * (UstarR - UR);
}

} // namespace cfd
