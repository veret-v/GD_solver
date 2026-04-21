#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "bc.hpp"
#include "riemann.hpp"
#include "mpi_parallel.hpp"

namespace cfd {

enum class SpatialScheme {
    Godunov,
    Kolgan
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

inline double minmod(double a, double b)
{
    if(a * b <= 0.0) return 0.0;
    return std::copysign(std::min(std::abs(a), std::abs(b)), a);
}

inline std::vector<std::array<Vec4, 2>> compute_gradients(const Mesh &m)
{
    std::vector<std::array<Vec4, 2>> grad(m.cells.size());
    for(int ci = 0; ci < m.nc(); ++ci) {
        const Cell &c = m.cells[ci];
        double a00 = 0.0;
        double a01 = 0.0;
        double a11 = 0.0;
        Vec4 b0 = make_vec4();
        Vec4 b1 = make_vec4();

        for(int fid : c.face_ids) {
            const int cj = other_cell(m.faces[fid], ci);
            if(cj < 0) continue;

            const double dx = m.cells[cj].cx - c.cx;
            const double dy = m.cells[cj].cy - c.cy;
            const double w = 1.0 / (dx * dx + dy * dy + 1e-30);
            a00 += w * dx * dx;
            a01 += w * dx * dy;
            a11 += w * dy * dy;

            for(int k = 0; k < 4; ++k) {
                const double du = m.cells[cj].U[k] - c.U[k];
                b0[k] += w * dx * du;
                b1[k] += w * dy * du;
            }
        }

        const double det = a00 * a11 - a01 * a01;
        if(std::abs(det) < 1e-30) continue;

        for(int k = 0; k < 4; ++k) {
            grad[ci][0][k] = (a11 * b0[k] - a01 * b1[k]) / det;
            grad[ci][1][k] = (a00 * b1[k] - a01 * b0[k]) / det;
        }
    }

    return grad;
}

inline Vec4 kolgan_state(const Mesh &m,
                         const std::vector<std::array<Vec4, 2>> &grad,
                         int cell_id,
                         int neighbor_id,
                         double x,
                         double y,
                         double g = FLOW.gamma)
{
    const Cell &c = m.cells[cell_id];
    const double dx = x - c.cx;
    const double dy = y - c.cy;
    Vec4 U = c.U;

    for(int k = 0; k < 4; ++k) {
        double slope = grad[cell_id][0][k] * dx + grad[cell_id][1][k] * dy;
        if(neighbor_id >= 0) {
            slope = minmod(slope, m.cells[neighbor_id].U[k] - c.U[k]);
        }
        U[k] += slope;
    }

    enforce_physical(U, g);
    return U;
}

// CHECK: UNSTRUCT_SCHEMES
inline void residuals_kolgan(Mesh &m, double g = FLOW.gamma)
{
    m.zero_res();
    const auto grad = compute_gradients(m);

    for(const Face &f : m.faces) {
        Vec4 UL = kolgan_state(m, grad, f.left, f.right, f.mx, f.my, g);
        Vec4 UR = f.is_boundary()
            ? ghost_state(UL, f.nx, f.ny, f.bc, g)
            : m.cells[f.right].U;

        if(f.right >= 0 && f.right < m.n_owned && f.bc == Face::BC::Interior) {
            UR = kolgan_state(m, grad, f.right, f.left, f.mx, f.my, g);
        }

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
    case SpatialScheme::Kolgan:
        residuals_kolgan(m, g);
        break;
    }
}

inline void update_cells(Mesh &m, double dt, double g = FLOW.gamma)
{
    for(int i = 0; i < m.n_owned; ++i) { // ВАЖНО: только n_owned!
        Cell &c = m.cells[i];
        c.U += (dt / c.vol) * c.res;
        enforce_physical(c.U, g);
    }
}

inline void step_euler(Mesh &m, Halo &halo, double dt, SpatialScheme scheme, double g = FLOW.gamma)
{
    exchange_halo(m, halo, MPI_COMM_WORLD); // Обмен перед расчетом
    compute_residuals(m, scheme, g);
    update_cells(m, dt, g);
}

inline void step_rk2(Mesh &m, double dt, SpatialScheme scheme, Halo& halo, MPI_Comm comm, double g = FLOW.gamma)
{
    std::vector<Vec4> U0(m.cells.size());
    for(size_t i=0; i<m.cells.size(); ++i) U0[i] = m.cells[i].U;

    // Стадия 1
    compute_residuals(m, scheme, g);
    update_cells(m, dt, g);
    
    // КРИТИЧЕСКИЙ ОБМЕН ТУТ
    exchange_halo(m, halo, comm); 

    // Стадия 2
    compute_residuals(m, scheme, g);
    for(int i=0; i < m.n_owned; ++i) {
        m.cells[i].U = 0.5 * (U0[i] + m.cells[i].U + (dt / m.cells[i].vol) * m.cells[i].res);
    }
    // После завершения шага тоже нужен обмен для следующего шага
    exchange_halo(m, halo, comm);
}

inline void step(Mesh &m, Halo &halo, double dt, SpatialScheme scheme, TimeScheme time_scheme, double g = FLOW.gamma)
{
    // Пробрасываем Halo
    switch(time_scheme) {
    case TimeScheme::Euler: step_euler(m, halo, dt, scheme, g); break;
    case TimeScheme::RK2:   step_rk2(m, dt, scheme, halo, MPI_COMM_WORLD, g); break;
    }
}

inline double compute_dt(const Mesh &m, double CFL, double g = FLOW.gamma)
{
    double dt = 1e100;
    for(int ci = 0; ci < m.n_owned; ++ci) { // Только n_owned!
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
    
    // Синхронизируем минимальный dt по всем процессам
    double global_dt;
    MPI_Allreduce(&dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global_dt;
}

inline Aero compute_aero(const Mesh& m, double rho_inf, double V_inf, double D_ref, double g=1.4) {
    double Fx=0, Fy=0, q=0.5*rho_inf*V_inf*V_inf;
    for (const auto& f : m.faces) {
        if (f.bc != Face::BC::Wall) continue;
        if (f.left >= m.n_owned) continue; // Учитываем стенку только на "родном" процессе

        const Vec4& U = m.cells[f.left].U;
        double rho=U[0], u=U[1]/rho, v=U[2]/rho;
        double p=(g-1)*(U[3]-0.5*rho*(u*u+v*v));
        Fx -= p * f.nx * f.length;
        Fy -= p * f.ny * f.length;
    }
    
    // Собираем силы со всех процессов
    double local_F[2] = {Fx, Fy}, global_F[2];
    MPI_Allreduce(local_F, global_F, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return {global_F[0]/(q*D_ref), global_F[1]/(q*D_ref)};
}

} // namespace cfd
