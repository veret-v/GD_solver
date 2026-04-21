#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

// Интерполяция для u-грани
double interpolate_u_to_face(const std::vector<std::vector<double>>& u,
                             const std::vector<std::vector<double>>& v,
                             int i, int j, char face);

double interpolate_v_to_face(const std::vector<std::vector<double>>& u,
                             const std::vector<std::vector<double>>& v,
                             int i, int j, char face);

// ========== СТАРЫЕ СИГНАТУРЫ (для обратной совместимости) ==========
// CHECK: PREDICTOR_U
void compute_momentum_coefficients_u(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& p,
    double dx, double dy, double nu, double rho,
    double dt, bool steady,
    std::vector<std::vector<double>>& aP,
    std::vector<std::vector<double>>& aE,
    std::vector<std::vector<double>>& aW,
    std::vector<std::vector<double>>& aN,
    std::vector<std::vector<double>>& aS,
    std::vector<std::vector<double>>& H,
    double U_lid
);
// CHECK: PREDICTOR_V
void compute_momentum_coefficients_v(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& p,
    double dx, double dy, double nu, double rho,
    double dt, bool steady,
    std::vector<std::vector<double>>& aP,
    std::vector<std::vector<double>>& aE,
    std::vector<std::vector<double>>& aW,
    std::vector<std::vector<double>>& aN,
    std::vector<std::vector<double>>& aS,
    std::vector<std::vector<double>>& H
);

void solve_pressure_correction(
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy,
    std::vector<std::vector<double>>& p_corr
);

void correct_velocity_pressure(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& p,
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& p_corr,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy, double alpha_p
);

// ========== НОВЫЕ СИГНАТУРЫ (с периодичностью и маской жидкости) ==========
void compute_momentum_coefficients_u(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& p,
    double dx, double dy, double nu, double rho,
    double dt, bool steady,
    std::vector<std::vector<double>>& aP,
    std::vector<std::vector<double>>& aE,
    std::vector<std::vector<double>>& aW,
    std::vector<std::vector<double>>& aN,
    std::vector<std::vector<double>>& aS,
    std::vector<std::vector<double>>& H,
    double U_lid,
    bool periodic_x,
    bool periodic_y,
    const std::vector<std::vector<int>>& fluid_mask
);

void compute_momentum_coefficients_v(
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& p,
    double dx, double dy, double nu, double rho,
    double dt, bool steady,
    std::vector<std::vector<double>>& aP,
    std::vector<std::vector<double>>& aE,
    std::vector<std::vector<double>>& aW,
    std::vector<std::vector<double>>& aN,
    std::vector<std::vector<double>>& aS,
    std::vector<std::vector<double>>& H,
    bool periodic_x,
    bool periodic_y,
    const std::vector<std::vector<int>>& fluid_mask
);

void solve_pressure_correction(
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy,
    std::vector<std::vector<double>>& p_corr,
    bool periodic_x,
    bool periodic_y,
    const std::vector<std::vector<int>>& fluid_mask
);

void correct_velocity_pressure(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& p,
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& p_corr,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy, double alpha_p,
    const std::vector<std::vector<int>>& fluid_mask
);

#endif // SOLVER_H