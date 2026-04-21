#include "solver.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

// ----------------------------------------------------------------------
// Вспомогательные функции интерполяции (без изменений)
// ----------------------------------------------------------------------
double interpolate_u_to_face(const std::vector<std::vector<double>>& u,
                             const std::vector<std::vector<double>>& v,
                             int i, int j, char face) {
    if (face == 'e') {
        return 0.5 * (u[i][j] + u[i+1][j]);
    } else if (face == 'w') {
        return 0.5 * (u[i][j] + u[i-1][j]);
    } else if (face == 'n') {
        return 0.5 * (v[i-1][j+1] + v[i][j+1]);
    } else if (face == 's') {
        return 0.5 * (v[i-1][j] + v[i][j]);
    }
    return 0.0;
}

double interpolate_v_to_face(const std::vector<std::vector<double>>& u,
                             const std::vector<std::vector<double>>& v,
                             int i, int j, char face) {
    if (face == 'e') {
        return 0.5 * (u[i+1][j-1] + u[i+1][j]);
    } else if (face == 'w') {
        return 0.5 * (u[i][j-1] + u[i][j]);
    } else if (face == 'n') {
        return 0.5 * (v[i][j] + v[i][j+1]);
    } else if (face == 's') {
        return 0.5 * (v[i][j] + v[i][j-1]);
    }
    return 0.0;
}

// ----------------------------------------------------------------------
// Функция вычисления коэффициентов для u-импульса (staggered, u на гранях x)
// ----------------------------------------------------------------------
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
    const std::vector<std::vector<int>>& fluid_mask)
{
    // Размеры сетки: u имеет размер (Nx+1) x Ny, где Nx - число ячеек по x
    int Nx = static_cast<int>(u.size()) - 1;
    int Ny = static_cast<int>(u[0].size());
    bool has_mask = !fluid_mask.empty();

    // Инициализация массивов коэффициентов (размер как у u)
    aP.assign(Nx + 1, std::vector<double>(Ny, 0.0));
    aE.assign(Nx + 1, std::vector<double>(Ny, 0.0));
    aW.assign(Nx + 1, std::vector<double>(Ny, 0.0));
    aN.assign(Nx + 1, std::vector<double>(Ny, 0.0));
    aS.assign(Nx + 1, std::vector<double>(Ny, 0.0));
    H.assign(Nx + 1, std::vector<double>(Ny, 0.0));

    double Dx = nu * dy / dx;
    double Dy = nu * dx / dy;

    for (int i = 1; i < Nx; ++i) {          // только внутренние u-грани
        for (int j = 0; j < Ny; ++j) {
            // --- Проверка маски (если обе соседние ячейки давления - стенка, пропускаем) ---
            if (has_mask) {
                bool left_fluid  = (i-1 >= 0 && j < Ny) ? (fluid_mask[i-1][j] == 1) : true;
                bool right_fluid = (i < Nx && j < Ny)   ? (fluid_mask[i][j] == 1)   : true;
                if (!left_fluid && !right_fluid) {
                    aP[i][j] = 1.0;
                    aE[i][j] = aW[i][j] = aN[i][j] = aS[i][j] = 0.0;
                    H[i][j] = 0.0;
                    continue;
                }
            }

            // --- Вычисление скоростей на гранях контрольного объёма (интерполяция) ---
            double ue = 0.5 * (u[i][j] + u[i+1][j]);
            double uw = 0.5 * (u[i][j] + u[i-1][j]);

            // v на северной и южной гранях (индексы v: i-1..i, j..j+1)
            // Индексы для v: v[iv][jv], где iv = i-1 или i, jv = j или j+1
            double vn, vs;
            if (periodic_y) {
                // Северная грань: jv = j+1, но если j == Ny-1, то периодически jv = 0
                int j_north = (j + 1) % Ny;
                vn = 0.5 * (v[i-1][j_north+1] + v[i][j_north+1]); // v имеет размер Nx x (Ny+1)
                vs = 0.5 * (v[i-1][j] + v[i][j]);
            } else {
                vn = 0.5 * (v[i-1][j+1] + v[i][j+1]);
                vs = 0.5 * (v[i-1][j]   + v[i][j]);
            }

            // Массовые потоки
            double Fe = rho * ue * dy;
            double Fw = rho * uw * dy;
            double Fn = rho * vn * dx;
            double Fs = rho * vs * dx;

            // Конвективно-диффузионные коэффициенты (схема upwind)
            aE[i][j] = std::max(-Fe, 0.0) + Dx;
            aW[i][j] = std::max( Fw, 0.0) + Dx;
            aN[i][j] = std::max(-Fn, 0.0) + Dy;
            aS[i][j] = std::max( Fs, 0.0) + Dy;

            // Граничные условия на стенках (если не периодические)
            if (j == 0 && !periodic_y) {
                // Нижняя стенка: расстояние dy/2
                aS[i][j] = std::max(Fs, 0.0) + nu * dx / (0.5 * dy);
            }
            if (j == Ny - 1 && !periodic_y) {
                // Верхняя стенка (движущаяся крышка или неподвижная)
                aN[i][j] = std::max(-Fn, 0.0) + nu * dx / (0.5 * dy);
            }

            // Диагональный коэффициент
            double ap_conv_diff = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
            aP[i][j] = ap_conv_diff;
            if (!steady) {
                aP[i][j] += rho * dx * dy / dt;
            }
            if (aP[i][j] < 1e-12) aP[i][j] = 1.0;  // защита от нуля

            // --- Значения u в соседних узлах для вычисления H ---
            double uE, uW, uN, uS;

            // Восточный сосед (i+1)
            if (periodic_x) {
                uE = (i < Nx - 1) ? u[i+1][j] : u[0][j];  // u[0] периодически равен u[Nx]
            } else {
                uE = (i < Nx - 1) ? u[i+1][j] : 0.0;
            }
            // Западный сосед (i-1)
            if (periodic_x) {
                uW = (i > 1) ? u[i-1][j] : u[Nx][j];      // u[Nx] периодически равен u[0]
            } else {
                uW = (i > 1) ? u[i-1][j] : 0.0;
            }
            // Северный сосед (j+1)
            if (periodic_y) {
                uN = (j < Ny - 1) ? u[i][j+1] : u[i][0];
            } else {
                uN = (j < Ny - 1) ? u[i][j+1] : 0.0;
            }
            // Южный сосед (j-1)
            if (periodic_y) {
                uS = (j > 0) ? u[i][j-1] : u[i][Ny-1];
            } else {
                uS = (j > 0) ? u[i][j-1] : 0.0;
            }

            // Ghost-значения для пристеночных ячеек (только для неподвижных стенок)
            if (!periodic_y) {
                if (j == 0) {
                    uS = -u[i][j];          // no-slip
                }
                if (j == Ny - 1) {
                    uN = 2.0 * U_lid - u[i][j];  // движущаяся крышка (если U_lid=0, то no-slip)
                }
            }

            // Сборка H
            H[i][j] = aE[i][j] * uE + aW[i][j] * uW + aN[i][j] * uN + aS[i][j] * uS;
            if (!steady) {
                H[i][j] += (rho * dx * dy / dt) * u[i][j];
            }
        }
    }
}

// ----------------------------------------------------------------------
// Функция вычисления коэффициентов для v-импульса (staggered, v на гранях y)
// ----------------------------------------------------------------------
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
    const std::vector<std::vector<int>>& fluid_mask)
{
    // Размеры: v имеет размер Nx x (Ny+1), где Ny - число ячеек по y
    int Nx = static_cast<int>(v.size());
    int Ny = static_cast<int>(v[0].size()) - 1;
    bool has_mask = !fluid_mask.empty();

    aP.assign(Nx, std::vector<double>(Ny + 1, 0.0));
    aE.assign(Nx, std::vector<double>(Ny + 1, 0.0));
    aW.assign(Nx, std::vector<double>(Ny + 1, 0.0));
    aN.assign(Nx, std::vector<double>(Ny + 1, 0.0));
    aS.assign(Nx, std::vector<double>(Ny + 1, 0.0));
    H.assign(Nx, std::vector<double>(Ny + 1, 0.0));

    double Dx = nu * dy / dx;
    double Dy = nu * dx / dy;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {  // только внутренние v-грани
            // Проверка маски
            if (has_mask) {
                bool bottom_fluid = (j-1 >= 0) ? (fluid_mask[i][j-1] == 1) : true;
                bool top_fluid    = (j < Ny)   ? (fluid_mask[i][j] == 1)   : true;
                if (!bottom_fluid && !top_fluid) {
                    aP[i][j] = 1.0;
                    aE[i][j] = aW[i][j] = aN[i][j] = aS[i][j] = 0.0;
                    H[i][j] = 0.0;
                    continue;
                }
            }

            // Скорости на гранях контрольного объёма v
            double ue, uw;
            if (periodic_x) {
                // Восточная грань: индекс i+1 (для u)
                int ie = (i < Nx - 1) ? i+1 : 0;
                ue = 0.5 * (u[ie][j-1] + u[ie][j]);
                // Западная грань: индекс i
                uw = 0.5 * (u[i][j-1] + u[i][j]);
            } else {
                ue = (i < Nx - 1) ? 0.5 * (u[i+1][j-1] + u[i+1][j]) : 0.0;
                uw = (i > 0)      ? 0.5 * (u[i][j-1]   + u[i][j])   : 0.0;
            }

            double vn = 0.5 * (v[i][j] + v[i][j+1]);
            double vs = 0.5 * (v[i][j] + v[i][j-1]);

            double Fe = rho * ue * dy;
            double Fw = rho * uw * dy;
            double Fn = rho * vn * dx;
            double Fs = rho * vs * dx;

            aE[i][j] = std::max(-Fe, 0.0) + Dx;
            aW[i][j] = std::max( Fw, 0.0) + Dx;
            aN[i][j] = std::max(-Fn, 0.0) + Dy;
            aS[i][j] = std::max( Fs, 0.0) + Dy;

            // Граничные условия на стенках
            if (i == 0 && !periodic_x) {
                aW[i][j] = std::max(Fw, 0.0) + nu * dy / (0.5 * dx);
            }
            if (i == Nx - 1 && !periodic_x) {
                aE[i][j] = std::max(-Fe, 0.0) + nu * dy / (0.5 * dx);
            }
            if (j == 1 && !periodic_y) {
                aS[i][j] = std::max(Fs, 0.0) + nu * dx / (0.5 * dy);
            }
            if (j == Ny - 1 && !periodic_y) {
                aN[i][j] = std::max(-Fn, 0.0) + nu * dx / (0.5 * dy);
            }

            double ap_conv_diff = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
            aP[i][j] = ap_conv_diff;
            if (!steady) {
                aP[i][j] += rho * dx * dy / dt;
            }
            if (aP[i][j] < 1e-12) aP[i][j] = 1.0;

            // Соседние значения v
            double vE, vW, vN, vS;

            if (periodic_x) {
                vE = (i < Nx - 1) ? v[i+1][j] : v[0][j];
                vW = (i > 0)      ? v[i-1][j] : v[Nx-1][j];
            } else {
                vE = (i < Nx - 1) ? v[i+1][j] : 0.0;
                vW = (i > 0)      ? v[i-1][j] : 0.0;
            }

            if (periodic_y) {
                vN = (j < Ny - 1) ? v[i][j+1] : v[i][0];    // v[i][0] периодически равен v[i][Ny]
                vS = (j > 1)      ? v[i][j-1] : v[i][Ny];   // v[i][Ny] периодически равен v[i][0]
            } else {
                vN = (j < Ny - 1) ? v[i][j+1] : 0.0;
                vS = (j > 1)      ? v[i][j-1] : 0.0;
            }

            // Ghost-значения для неподвижных стенок
            if (!periodic_x) {
                if (i == 0)        vW = -v[i][j];
                if (i == Nx - 1)   vE = -v[i][j];
            }
            if (!periodic_y) {
                if (j == 1)        vS = -v[i][j];
                if (j == Ny - 1)   vN = -v[i][j];
            }

            H[i][j] = aE[i][j] * vE + aW[i][j] * vW + aN[i][j] * vN + aS[i][j] * vS;
            if (!steady) {
                H[i][j] += (rho * dx * dy / dt) * v[i][j];
            }
        }
    }
}

// ----------------------------------------------------------------------
// НОВАЯ РЕАЛИЗАЦИЯ: решение уравнения Пуассона для p' (с периодичностью и маской)
// ----------------------------------------------------------------------
void solve_pressure_correction(
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy,
    std::vector<std::vector<double>>& p_corr,
    bool periodic_x,
    bool periodic_y,
    const std::vector<std::vector<int>>& fluid_mask)
{
    int Nx = p_corr.size();
    int Ny = p_corr[0].size();
    bool has_mask = !fluid_mask.empty();

    std::vector<std::vector<double>> aE(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> aW(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> aN(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> aS(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> b(Nx, std::vector<double>(Ny, 0.0));

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if (has_mask && fluid_mask[i][j] == 0) {
                aE[i][j] = aW[i][j] = aN[i][j] = aS[i][j] = 0.0;
                b[i][j] = 0.0;
                continue;
            }

            // Коэффициенты для соседних ячеек
            if (i < Nx-1) {
                aE[i][j] = (dy * dy) / aP_u[i+1][j];
            } else if (periodic_x) {
                aE[i][j] = (dy * dy) / aP_u[0][j];
            }

            if (i > 0) {
                aW[i][j] = (dy * dy) / aP_u[i][j];
            } else if (periodic_x) {
                aW[i][j] = (dy * dy) / aP_u[Nx-1][j];
            }

            if (j < Ny-1) {
                aN[i][j] = (dx * dx) / aP_v[i][j+1];
            } else if (periodic_y) {
                aN[i][j] = (dx * dx) / aP_v[i][0];
            }

            if (j > 0) {
                aS[i][j] = (dx * dx) / aP_v[i][j];
            } else if (periodic_y) {
                aS[i][j] = (dx * dx) / aP_v[i][Ny-1];
            }

            double ue = u_star[i+1][j];
            double uw = u_star[i][j];
            double vn = v_star[i][j+1];
            double vs = v_star[i][j];
            double div = (ue - uw) * dy + (vn - vs) * dx;
            b[i][j] = -div;
        }
    }

    bool fix_pressure = !(periodic_x && periodic_y);

    double omega = 1.5;
    int max_iter = 5000;
    double tol = 1e-8;

    for (int iter = 0; iter < max_iter; ++iter) {
        double res_norm = 0.0;
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                if (has_mask && fluid_mask[i][j] == 0) {
                    p_corr[i][j] = 0.0;
                    continue;
                }
                if (fix_pressure && i == 0 && j == 0) {
                    p_corr[i][j] = 0.0;
                    continue;
                }

                double ae = aE[i][j];
                double aw = aW[i][j];
                double an = aN[i][j];
                double as = aS[i][j];
                double ap = ae + aw + an + as;
                if (ap < 1e-12) ap = 1e12;

                double p_e = (i < Nx-1) ? p_corr[i+1][j] : (periodic_x ? p_corr[0][j] : 0.0);
                double p_w = (i > 0)    ? p_corr[i-1][j] : (periodic_x ? p_corr[Nx-1][j] : 0.0);
                double p_n = (j < Ny-1) ? p_corr[i][j+1] : (periodic_y ? p_corr[i][0] : 0.0);
                double p_s = (j > 0)    ? p_corr[i][j-1] : (periodic_y ? p_corr[i][Ny-1] : 0.0);

                double p_new = (ae * p_e + aw * p_w + an * p_n + as * p_s + b[i][j]) / ap;
                p_corr[i][j] = p_corr[i][j] + omega * (p_new - p_corr[i][j]);

                double res = ap * p_corr[i][j] - (ae * p_e + aw * p_w + an * p_n + as * p_s + b[i][j]);
                res_norm += res * res;
            }
        }
        if (sqrt(res_norm / (Nx*Ny)) < tol) break;
        
    }

}

// ----------------------------------------------------------------------
// НОВАЯ РЕАЛИЗАЦИЯ: коррекция скорости и давления (с маской)
// ----------------------------------------------------------------------
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
    const std::vector<std::vector<int>>& fluid_mask)
{
    int Nx_u = u.size() - 1;
    int Ny_u = u[0].size();
    int Nx_v = v.size();
    int Ny_v = v[0].size() - 1;
    int Nx_p = p.size();
    int Ny_p = p[0].size();
    bool has_mask = !fluid_mask.empty();

    // Давление
    for (int i = 0; i < Nx_p; ++i) {
        for (int j = 0; j < Ny_p; ++j) {
            if (has_mask && fluid_mask[i][j] == 0) continue;
            p[i][j] += alpha_p * p_corr[i][j];
        }
    }

    // Скорость u
    for (int i = 1; i < Nx_u; ++i) {
        for (int j = 0; j < Ny_u; ++j) {
            if (has_mask) {
                bool left_fluid  = (i-1 >= 0 && j < Ny_u) ? (fluid_mask[i-1][j] == 1) : true;
                bool right_fluid = (i < Nx_p && j < Ny_u) ? (fluid_mask[i][j] == 1)   : true;
                if (!left_fluid && !right_fluid) {
                    u[i][j] = 0.0;
                    continue;
                }
            }
            double dp = p_corr[i][j] - p_corr[i-1][j];
            u[i][j] = u_star[i][j] - (dy / aP_u[i][j]) * dp;
        }
    }

    // Скорость v
    for (int i = 0; i < Nx_v; ++i) {
        for (int j = 1; j < Ny_v; ++j) {
            if (has_mask) {
                bool bottom_fluid = (j-1 >= 0) ? (fluid_mask[i][j-1] == 1) : true;
                bool top_fluid    = (j < Ny_p) ? (fluid_mask[i][j] == 1)   : true;
                if (!bottom_fluid && !top_fluid) {
                    v[i][j] = 0.0;
                    continue;
                }
            }
            double dp = p_corr[i][j] - p_corr[i][j-1];
            v[i][j] = v_star[i][j] - (dx / aP_v[i][j]) * dp;
        }
    }
}

// ----------------------------------------------------------------------
// СТАРЫЕ СИГНАТУРЫ (обёртки для обратной совместимости)
// ----------------------------------------------------------------------
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
    double U_lid)
{
    std::vector<std::vector<int>> empty_mask;
    compute_momentum_coefficients_u(u, v, p, dx, dy, nu, rho, dt, steady,
                                    aP, aE, aW, aN, aS, H, U_lid,
                                    false, false, empty_mask);
}

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
    std::vector<std::vector<double>>& H)
{
    std::vector<std::vector<int>> empty_mask;
    compute_momentum_coefficients_v(u, v, p, dx, dy, nu, rho, dt, steady,
                                    aP, aE, aW, aN, aS, H,
                                    false, false, empty_mask);
}

void solve_pressure_correction(
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy,
    std::vector<std::vector<double>>& p_corr)
{
    std::vector<std::vector<int>> empty_mask;
    solve_pressure_correction(u_star, v_star, aP_u, aP_v, dx, dy, p_corr,
                              false, false, empty_mask);
}

void correct_velocity_pressure(
    std::vector<std::vector<double>>& u,
    std::vector<std::vector<double>>& v,
    std::vector<std::vector<double>>& p,
    const std::vector<std::vector<double>>& u_star,
    const std::vector<std::vector<double>>& v_star,
    const std::vector<std::vector<double>>& p_corr,
    const std::vector<std::vector<double>>& aP_u,
    const std::vector<std::vector<double>>& aP_v,
    double dx, double dy, double alpha_p)
{
    std::vector<std::vector<int>> empty_mask;
    correct_velocity_pressure(u, v, p, u_star, v_star, p_corr, aP_u, aP_v,
                              dx, dy, alpha_p, empty_mask);
}