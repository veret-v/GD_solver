// main.cpp — последовательная реализация SIMPLE / PISO / PIMPLE
// Поддержка тестов: cavity, taylor_green, backward_step
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <filesystem>
#include "./parser.h"
#include "./solver.h"

// Структуры параметров для каждого метода
struct SIMPLEParams {
    double alpha_u;
    double alpha_p;
    int max_iter;
    double tol;
};

struct PISOParams {
    double dt;
    double t_end;
    int n_correctors;
};

struct PIMPLEParams {
    double dt;
    double t_end;
    int n_outer;
    int n_correctors;
    double alpha_u;
    double alpha_p;
};

// Структура для описания тестового случая
struct TestCase {
    std::string name;
    double Lx, Ly;
    double U_lid;               // для каверны
    double U_inlet;             // для ступеньки (скорость на входе)
    double nu;                  // кинематическая вязкость
    double Re;                  // число Рейнольдса
    double t_end;               // конечное время (для нестационара)
    double dt;                  // шаг по времени
    int Nx, Ny;
    bool is_steady;             // стационарная задача или нет
    std::vector<std::vector<int>> fluid_mask; // 1 – жидкость, 0 – стенка (для ступеньки)
};

// ----------------------------------------------------------------------
// Вспомогательные функции инициализации и граничных условий
// ----------------------------------------------------------------------

// Инициализация полей для вихря Тейлора–Грина
// CHECK: TG_INIT
void init_taylor_green(std::vector<std::vector<double>>& u,
                       std::vector<std::vector<double>>& v,
                       std::vector<std::vector<double>>& p,
                       double Lx, double Ly, int Nx, int Ny, double nu) {
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    
    // u на гранях (Nx+1) x Ny
    for (int i = 0; i <= Nx; ++i) {
        double x = i * dx;               // u-грани: i=0 => x=0; i=Nx => x=Lx
        for (int j = 0; j < Ny; ++j) {
            double y = (j + 0.5) * dy;   // центр ячейки по y
            u[i][j] = std::sin(x) * std::cos(y);
        }
    }
    
    // v на гранях Nx x (Ny+1)
    for (int i = 0; i < Nx; ++i) {
        double x = (i + 0.5) * dx;       // центр ячейки по x
        for (int j = 0; j <= Ny; ++j) {
            double y = j * dy;           // v-грани: j=0 => y=0; j=Ny => y=Ly
            v[i][j] = -std::cos(x) * std::sin(y);
        }
    }
    
    // Давление в центрах ячеек
    for (int i = 0; i < Nx; ++i) {
        double x = (i + 0.5) * dx;
        for (int j = 0; j < Ny; ++j) {
            double y = (j + 0.5) * dy;
            p[i][j] = -0.25 * (std::cos(2*x) + std::cos(2*y));
        }
    }
}

// Создание маски жидкости для обратной ступеньки
std::vector<std::vector<int>> create_backward_step_mask(int Nx, int Ny, double Lx, double Ly) {
    std::vector<std::vector<int>> mask(Nx, std::vector<int>(Ny, 1)); // по умолчанию жидкость
    double step_height = Ly / 2.0;
    double step_length = Lx / 6.0;   // длина ступеньки 1/6 области (можно настроить)
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    for (int i = 0; i < Nx; ++i) {
        double x = (i + 0.5) * dx;
        for (int j = 0; j < Ny; ++j) {
            double y = (j + 0.5) * dy;
            if (x < step_length && y < step_height) {
                mask[i][j] = 0; // стенка
            }
        }
    }
    return mask;
}

// Граничные условия для каверны (стандартные)
void apply_boundary_conditions_cavity(std::vector<std::vector<double>>& u,
                                      std::vector<std::vector<double>>& v,
                                      int Nx, int Ny, double U_lid) {
    // u: левая и правая стенки
    for (int j = 0; j < Ny; ++j) u[0][j] = 0.0;
    for (int j = 0; j < Ny; ++j) u[Nx][j] = 0.0;
    // верхняя крышка
    for (int i = 0; i <= Nx; ++i) u[i][Ny-1] = U_lid;
    // v: нижняя и верхняя стенки
    for (int i = 0; i < Nx; ++i) v[i][0] = 0.0;
    for (int i = 0; i < Nx; ++i) v[i][Ny] = 0.0;
    // v: боковые стенки
    for (int j = 1; j < Ny; ++j) {
        v[0][j] = 0.0;
        v[Nx-1][j] = 0.0;
    }
}


// CHECK: TG_BC
// Новые граничные условия Дирихле: точное решение на всех границах
void apply_boundary_conditions_taylor_green(std::vector<std::vector<double>>& u,
                                            std::vector<std::vector<double>>& v,
                                            int Nx, int Ny,
                                            double Lx, double Ly,
                                            double nu,
                                            double time)   // время для нестационара
{
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    
    // Аналитическое решение вихря Тейлора–Грина:
    // u(x,y,t) =  sin(x) * cos(y) * exp(-2 nu t)
    // v(x,y,t) = -cos(x) * sin(y) * exp(-2 nu t)
    double factor = std::exp(-2.0 * nu * time);
    
    // ---------- u-границы ----------
    // Левая граница i = 0
    for (int j = 0; j < Ny; ++j) {
        double y = (j + 0.5) * dy;
        double x = 0.0;                     // u-грань на левой границе
        u[0][j] = std::sin(x) * std::cos(y) * factor;
    }
    // Правая граница i = Nx
    for (int j = 0; j < Ny; ++j) {
        double y = (j + 0.5) * dy;
        double x = Lx;                      // u-грань на правой границе
        u[Nx][j] = std::sin(x) * std::cos(y) * factor;
    }
    // Верхняя и нижняя границы для u: значения в центрах ячеек (i, j=0 и j=Ny-1)
    // Эти грани не являются границами для u на разнесённой сетке, но для полноты:
    // (на самом деле u не хранится на y-границах, поэтому не трогаем)
    
    // ---------- v-границы ----------
    // Нижняя граница j = 0
    for (int i = 0; i < Nx; ++i) {
        double x = (i + 0.5) * dx;
        double y = 0.0;                     // v-грань на нижней границе
        v[i][0] = -std::cos(x) * std::sin(y) * factor;
    }
    // Верхняя граница j = Ny
    for (int i = 0; i < Nx; ++i) {
        double x = (i + 0.5) * dx;
        double y = Ly;                      // v-грань на верхней границе
        v[i][Ny] = -std::cos(x) * std::sin(y) * factor;
    }
    // Левая и правая границы для v: значения в центрах (i=0 и i=Nx-1)
    // (аналогично, не являются границами для v)
}


// Граничные условия для обратной ступеньки (входной профиль, выход, стенки)
// CHECK: STEP_BLOCK
void apply_boundary_conditions_step(std::vector<std::vector<double>>& u,
                                    std::vector<std::vector<double>>& v,
                                    int Nx, int Ny,
                                    double U_inlet,
                                    const std::vector<std::vector<int>>& fluid_mask,
                                    double Ly) {
    double dy = Ly / Ny;
    // Вход слева (i=0): параболический профиль u(0,y) = 4*U_inlet * (y/Ly)*(1 - y/Ly)
    // CHECK: STEP_INLET
    for (int j = 0; j < Ny; ++j) {
        double y = (j + 0.5) * dy;
        double u_in = 6.0 * U_inlet * (y - Ly / 2) * (Ly - y) / (Ly - Ly / 2);
        u[0][j] = u_in;
    }
    // Выход справа (i=Nx): нулевой градиент — экстраполяция
    for (int j = 0; j < Ny; ++j) {
        u[Nx][j] = u[Nx-1][j];
    }
    // Верхняя и нижняя стенки для v
    for (int i = 0; i < Nx; ++i) {
        v[i][0] = 0.0;
        v[i][Ny] = 0.0;
    }
    // Боковые стенки для v (кроме входа/выхода)
    for (int j = 1; j < Ny; ++j) {
        v[0][j] = 0.0;
        v[Nx-1][j] = 0.0;
    }
    // Если есть маска, зануляем скорость в твёрдых ячейках
    if (!fluid_mask.empty()) {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                if (fluid_mask[i][j] == 0) {
                    // u на правой и левой гранях ячейки
                    if (i+1 <= Nx) u[i+1][j] = 0.0;
                    if (i > 0) u[i][j] = 0.0;
                    // v на верхней и нижней гранях
                    if (j+1 <= Ny) v[i][j+1] = 0.0;
                    if (j > 0) v[i][j] = 0.0;
                }
            }
        }
    }
}

// ----------------------------------------------------------------------
// Основная функция для запуска одного теста
// ----------------------------------------------------------------------
void run_single_test(const TestCase& tcase,
                     int equation_type,
                     const SIMPLEParams& simp,
                     const PISOParams& piso,
                     const PIMPLEParams& pimple,
                     int step_out) {
    std::cout << "\n===== Running test: " << tcase.name << " =====\n";
    double dx = tcase.Lx / tcase.Nx;
    double dy = tcase.Ly / tcase.Ny;
    double rho = 1.0;
    double nu = tcase.nu;

    // Инициализация полей
    // CHECK: STAGGERED_GRID
    std::vector<std::vector<double>> u(tcase.Nx+1, std::vector<double>(tcase.Ny, 0.0));
    std::vector<std::vector<double>> v(tcase.Nx, std::vector<double>(tcase.Ny+1, 0.0));
    std::vector<std::vector<double>> p(tcase.Nx, std::vector<double>(tcase.Ny, 0.0));

    bool periodic_x = false;
    bool periodic_y = false;
    const auto& mask = tcase.fluid_mask;

    // Установка начальных условий
    if (tcase.name == "cavity") {
        // всё уже нули
    } else if (tcase.name == "taylor_green") {
        init_taylor_green(u, v, p, tcase.Lx, tcase.Ly, tcase.Nx, tcase.Ny, nu);
    } else if (tcase.name == "backward_step") {
        // нулевые начальные скорости
    }

    // Применение граничных условий (начальное)
    if (tcase.name == "cavity") {
        apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
    } else if (tcase.name == "taylor_green") {
        apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, 0.0);
    } else if (tcase.name == "backward_step") {
        apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);
    }

    // Создание директории для вывода
    std::string test_output_dir = "../../output/" + tcase.name + "/";
    std::filesystem::create_directories(test_output_dir);
    auto save_state_test = [&](int step, double time, const std::string& suffix = "") {
        std::ostringstream filename;
        filename << test_output_dir << "step_" << std::setw(5) << std::setfill('0') << step
                 << "_t_" << std::fixed << std::setprecision(4) << time << suffix << ".csv";
        std::ofstream fout(filename.str());
        fout << "x,y,u,v,p\n";
        for (int i = 0; i < tcase.Nx; ++i) {
            double x = (i + 0.5) * dx;
            for (int j = 0; j < tcase.Ny; ++j) {
                double y = (j + 0.5) * dy;
                double u_c = 0.5 * (u[i][j] + u[i+1][j]);
                double v_c = 0.5 * (v[i][j] + v[i][j+1]);
                fout << x << "," << y << "," << u_c << "," << v_c << "," << p[i][j] << "\n";
            }
        }
        fout.close();
    };

    save_state_test(0, 0.0, "_initial");

    auto start_time = std::chrono::high_resolution_clock::now();

    // ------------------------------------------------------------------
    // Расчётный блок
    // ------------------------------------------------------------------
    if (equation_type == 8) { // SIMPLE (steady)
        std::cout << "SIMPLE: max_iter=" << simp.max_iter << ", tol=" << simp.tol << "\n";
        // CHECK: SIMPLE_LOOP
        for (int iter = 1; iter <= simp.max_iter; ++iter) {
            auto u_old = u;
            auto v_old = v;

            // Коэффициенты и предиктор u*
            // CHECK: PREDICTOR_U
            std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
            compute_momentum_coefficients_u(u, v, p, dx, dy, nu, rho, 0.0, true,
                                aP_u, aE_u, aW_u, aN_u, aS_u, H_u,
                                tcase.U_lid, periodic_x, periodic_y, mask);
            // Граничные условия
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, 0.0);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            std::vector<std::vector<double>> u_star = u;
            for (int i = 1; i < tcase.Nx; ++i) {
                for (int j = 0; j < tcase.Ny; ++j) {
                    double H = H_u[i][j];
                    double dp = (p[i][j] - p[i-1][j]) * dy;
                    u_star[i][j] = (1.0 - simp.alpha_u) * u[i][j] +
                                    (simp.alpha_u / aP_u[i][j]) * (H - dp);
                }
            }

            // Коэффициенты и предиктор v*
            // CHECK: PREDICTOR_V
            std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
            compute_momentum_coefficients_v(u, v, p, dx, dy, nu, rho, 0.0, true,
                                aP_v, aE_v, aW_v, aN_v, aS_v, H_v,
                                periodic_x, periodic_y, mask);
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, 0.0);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            std::vector<std::vector<double>> v_star = v;
            for (int i = 0; i < tcase.Nx; ++i) {
                for (int j = 1; j < tcase.Ny; ++j) {
                    double H = H_v[i][j];
                    double dp = (p[i][j] - p[i][j-1]) * dx;
                    v_star[i][j] = (1.0 - simp.alpha_u) * v[i][j] +
                                    (simp.alpha_u / aP_v[i][j]) * (H - dp);
                }
            }

            // Уравнение Пуассона для p'
            // CHECK: POISSON
            std::vector<std::vector<double>> p_corr(tcase.Nx, std::vector<double>(tcase.Ny, 0.0));
            solve_pressure_correction(u_star, v_star, aP_u, aP_v, dx, dy, p_corr,
                          periodic_x, periodic_y, mask);

            // Коррекция
            // CHECK: CORRECTION
            correct_velocity_pressure(u, v, p, u_star, v_star, p_corr, aP_u, aP_v, dx, dy, simp.alpha_p, mask);
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, 0.0);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            // Проверка сходимости
            double max_du = 0.0;
            for (int i = 1; i < tcase.Nx; ++i)
                for (int j = 0; j < tcase.Ny; ++j)
                    max_du = std::max(max_du, std::abs(u[i][j] - u_old[i][j]));
            for (int i = 0; i < tcase.Nx; ++i)
                for (int j = 1; j < tcase.Ny; ++j)
                    max_du = std::max(max_du, std::abs(v[i][j] - v_old[i][j]));

            if (iter % 100 == 0 || max_du < simp.tol) {
                std::cout << "  Iter " << iter << ", max_du = " << max_du << "\n";
                if (iter % step_out == 0 || max_du < simp.tol)
                    save_state_test(iter, iter, "");
            }
            if (max_du < simp.tol) {
                std::cout << "  SIMPLE converged at iteration " << iter << "\n";
                break;
            }
        }
    }
    else if (equation_type == 9) { // PISO (unsteady)
        std::cout << "PISO: dt=" << piso.dt << ", t_end=" << piso.t_end
                  << ", n_correctors=" << piso.n_correctors << "\n";
        double t = 0.0;
        int step = 0;
        save_state_test(0, 0.0, "");
        // CHECK: PISO_LOOP
        while (t < piso.t_end) {
            auto u_old = u;
            auto v_old = v;
            auto p_old = p;

            // Предиктор u*
            // CHECK: PREDICTOR_U
            std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
            compute_momentum_coefficients_u(u_old, v_old, p_old, dx, dy, nu, rho,
                                            piso.dt, false,   // <-- нестационарный
                                            aP_u, aE_u, aW_u, aN_u, aS_u, H_u,
                                            tcase.U_lid, periodic_x, periodic_y, mask);
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            std::vector<std::vector<double>> u_star = u_old;
            for (int i = 1; i < tcase.Nx; ++i) {
                for (int j = 0; j < tcase.Ny; ++j) {
                    double H = H_u[i][j];
                    double dp = (p_old[i][j] - p_old[i-1][j]) * dy;
                    u_star[i][j] = (H - dp) / aP_u[i][j];
                }
            }

            // Предиктор v*
            // CHECK: PREDICTOR_V
            std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
            compute_momentum_coefficients_v(u_old, v_old, p_old, dx, dy, nu, rho,
                                            piso.dt, false,   // <-- нестационарный
                                            aP_v, aE_v, aW_v, aN_v, aS_v, H_v,
                                            periodic_x, periodic_y, mask);
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            std::vector<std::vector<double>> v_star = v_old;
            for (int i = 0; i < tcase.Nx; ++i) {
                for (int j = 1; j < tcase.Ny; ++j) {
                    double H = H_v[i][j];
                    double dp = (p_old[i][j] - p_old[i][j-1]) * dx;
                    v_star[i][j] = (H - dp) / aP_v[i][j];
                }
            }

            u = u_star;
            v = v_star;
            p = p_old;

            // Коррекции
            for (int corr = 0; corr < piso.n_correctors; ++corr) {
                std::vector<std::vector<double>> p_corr(tcase.Nx, std::vector<double>(tcase.Ny, 0.0));
                // CHECK: POISSON
                solve_pressure_correction(u, v, aP_u, aP_v, dx, dy, p_corr,
                          periodic_x, periodic_y, mask);
                // CHECK: CORRECTION
                correct_velocity_pressure(u, v, p, u, v, p_corr, aP_u, aP_v, dx, dy, 1.0, mask);
            }
            if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
            else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
            else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

            t += piso.dt;
            step++;
            if (step % step_out == 0) {
                std::cout << "  Step " << step << ", t = " << t << "\n";
                save_state_test(step, t, "");
            }
        }
    }
    else if (equation_type == 10) { // PIMPLE (unsteady)
        std::cout << "PIMPLE: dt=" << pimple.dt << ", t_end=" << pimple.t_end
                  << ", nOuter=" << pimple.n_outer << ", nCorr=" << pimple.n_correctors << "\n";
        double t = 0.0;
        int step = 0;
        save_state_test(0, 0.0, "");
        // CHECK: PIMPLE_LOOP
        while (t < pimple.t_end) {
            auto u_old = u;
            auto v_old = v;
            auto p_old = p;

            for (int outer = 0; outer < pimple.n_outer; ++outer) {
                double alpha_u_cur = (outer == pimple.n_outer - 1) ? 1.0 : pimple.alpha_u;
                double alpha_p_cur = (outer == pimple.n_outer - 1) ? 1.0 : pimple.alpha_p;

                // u
                // CHECK: PREDICTOR_U
                std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
                compute_momentum_coefficients_u(u, v, p, dx, dy, nu, rho,
                                                pimple.dt, false,   // <-- нестационарный
                                                aP_u, aE_u, aW_u, aN_u, aS_u, H_u,
                                                tcase.U_lid, periodic_x, periodic_y, mask);
                if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
                else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
                else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

                std::vector<std::vector<double>> u_star = u;
                for (int i = 1; i < tcase.Nx; ++i) {
                    for (int j = 0; j < tcase.Ny; ++j) {
                        double H = H_u[i][j];
                        double dp = (p[i][j] - p[i-1][j]) * dy;
                        u_star[i][j] = (1.0 - alpha_u_cur) * u[i][j] +
                                        (alpha_u_cur / aP_u[i][j]) * (H - dp);
                        aP_u[i][j] /= alpha_u_cur;
                    }
                }

                // v
                // CHECK: PREDICTOR_V
                std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
                compute_momentum_coefficients_v(u, v, p, dx, dy, nu, rho,
                                                pimple.dt, false,   // <-- нестационарный
                                                aP_v, aE_v, aW_v, aN_v, aS_v, H_v,
                                                periodic_x, periodic_y, mask);
                if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
                else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
                else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

                std::vector<std::vector<double>> v_star = v;
                for (int i = 0; i < tcase.Nx; ++i) {
                    for (int j = 1; j < tcase.Ny; ++j) {
                        double H = H_v[i][j];
                        double dp = (p[i][j] - p[i][j-1]) * dx;
                        v_star[i][j] = (1.0 - alpha_u_cur) * v[i][j] +
                                        (alpha_u_cur / aP_v[i][j]) * (H - dp);
                        aP_v[i][j] /= alpha_u_cur;
                    }
                }

                std::vector<std::vector<double>> u_corr = u_star;
                std::vector<std::vector<double>> v_corr = v_star;
                std::vector<std::vector<double>> p_corr_field = p;
                for (int corr = 0; corr < pimple.n_correctors; ++corr) {
                    std::vector<std::vector<double>> p_corr(tcase.Nx, std::vector<double>(tcase.Ny, 0.0));
                    // CHECK: POISSON
                    solve_pressure_correction(u_corr, v_corr, aP_u, aP_v, dx, dy, p_corr,
                          periodic_x, periodic_y, mask);
                    // CHECK: CORRECTION
                    correct_velocity_pressure(u_corr, v_corr, p_corr_field,
                                              u_corr, v_corr, p_corr,
                                              aP_u, aP_v, dx, dy, alpha_p_cur, mask);
                }
                if (tcase.name == "cavity") apply_boundary_conditions_cavity(u, v, tcase.Nx, tcase.Ny, tcase.U_lid);
                else if (tcase.name == "taylor_green") apply_boundary_conditions_taylor_green(u, v, tcase.Nx, tcase.Ny, tcase.Lx, tcase.Ly, nu, t);
                else if (tcase.name == "backward_step") apply_boundary_conditions_step(u, v, tcase.Nx, tcase.Ny, tcase.U_inlet, tcase.fluid_mask, tcase.Ly);

                u = u_corr;
                v = v_corr;
                p = p_corr_field;
            }

            t += pimple.dt;
            step++;
            if (step % step_out == 0) {
                std::cout << "  Step " << step << ", t = " << t << "\n";
                save_state_test(step, t, "");
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Test " << tcase.name << " completed in " << duration.count() / 1000.0 << " s.\n";
    save_state_test(9999, (equation_type == 8) ? simp.max_iter : tcase.t_end, "_final");
}

// ----------------------------------------------------------------------
// MAIN
// ----------------------------------------------------------------------
int main(int argc, char** argv) {
    // Чтение конфигурации
    std::string system_ini = "../../configs/input.ini";
    IniParser sys(system_ini);

    // Общие параметры
    std::string case_name = sys.getString("case_name");
    double x_min = sys.getDouble("x_min");
    double x_max = sys.getDouble("x_max");
    double y_min = sys.getDouble("y_min");
    double y_max = sys.getDouble("y_max");
    int Nx = sys.getInt("Nx");
    int Ny = sys.getInt("Ny");
    double tmax = sys.getDouble("tmax");
    double cfl = sys.getDouble("cfl");
    double Re = sys.getDouble("Re");
    double U_lid = sys.getDouble("U_lid");
    int equation_type = sys.getInt("equation_type"); // 8=SIMPLE, 9=PISO, 10=PIMPLE

    // Новые параметры для тестов
    std::string test_case_name = sys.getString("test_case");
    int run_all = sys.getInt("run_all_tests");
    double U_inlet = sys.getDouble("U_inlet");   // скорость на входе для ступеньки

    // Параметры методов
    SIMPLEParams simp;
    PISOParams piso;
    PIMPLEParams pimple;

    if (equation_type == 8) {
        simp.alpha_u = sys.getDouble("alpha_u");
        simp.alpha_p = sys.getDouble("alpha_p");
        simp.max_iter = sys.getInt("max_iter_simple");
        simp.tol = sys.getDouble("tol_simple");
    } else if (equation_type == 9) {
        piso.dt = sys.getDouble("dt");
        piso.t_end = tmax;
        piso.n_correctors = sys.getInt("n_correctors");
    } else if (equation_type == 10) {
        pimple.dt = sys.getDouble("dt");
        pimple.t_end = tmax;
        pimple.n_outer = sys.getInt("nOuterCorr");
        pimple.n_correctors = sys.getInt("nCorrectors");
        pimple.alpha_u = sys.getDouble("alpha_u");
        pimple.alpha_p = sys.getDouble("alpha_p");
    }

    int step_out = sys.getInt("step_out");

    // Формирование списка тестов
    std::vector<TestCase> tests;

    if (run_all) {
        // Тест 1: каверна
        TestCase cav;
        cav.name = "cavity";
        cav.Lx = x_max - x_min;
        cav.Ly = y_max - y_min;
        cav.Nx = Nx;
        cav.Ny = Ny;
        cav.Re = Re;
        cav.nu = 1.0 / Re;
        cav.U_lid = U_lid;
        cav.U_inlet = 0.0;
        cav.t_end = tmax;
        cav.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.0;
        cav.is_steady = (equation_type == 8);
        tests.push_back(cav);

        // Тест 2: вихрь Тейлора–Грина
        TestCase tg;
        tg.name = "taylor_green";
        tg.Lx = 2.0 * M_PI;
        tg.Ly = 2.0 * M_PI;
        tg.Nx = Nx;
        tg.Ny = Ny;
        tg.Re = Re;
        tg.nu = 1.0 / Re;
        tg.U_lid = 0.0;
        tg.U_inlet = 0.0;
        tg.t_end = tmax;
        tg.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.01;
        tg.is_steady = false;
        tests.push_back(tg);

        // Тест 3: обратная ступенька
        TestCase bs;
        bs.name = "backward_step";
        bs.Lx = 15.0;          // длина канала
        bs.Ly = 1.0;           // высота канала
        bs.Nx = Nx;
        bs.Ny = Ny;
        bs.Re = Re;
        bs.nu = 1.0 / Re;
        bs.U_lid = 0.0;
        bs.U_inlet = U_inlet;
        bs.t_end = tmax;
        bs.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.01;
        bs.is_steady = (equation_type == 8);
        bs.fluid_mask = create_backward_step_mask(bs.Nx, bs.Ny, bs.Lx, bs.Ly);
        tests.push_back(bs);
    } else {
        // Одиночный тест
        TestCase single;
        single.name = test_case_name;
        if (test_case_name == "cavity") {
            single.Lx = x_max - x_min;
            single.Ly = y_max - y_min;
            single.Nx = Nx;
            single.Ny = Ny;
            single.Re = Re;
            single.nu = 1.0 / Re;
            single.U_lid = U_lid;
            single.U_inlet = 0.0;
            single.t_end = tmax;
            single.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.0;
            single.is_steady = (equation_type == 8);
        } else if (test_case_name == "taylor_green") {
            single.Lx = 2.0 * M_PI;
            single.Ly = 2.0 * M_PI;
            single.Nx = Nx;
            single.Ny = Ny;
            single.Re = Re;
            single.nu = 1.0 / Re;
            single.U_lid = 0.0;
            single.U_inlet = 0.0;
            single.t_end = tmax;
            single.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.01;
            single.is_steady = false;
        } else if (test_case_name == "backward_step") {
            single.Lx = 15.0;
            single.Ly = 1.0;
            single.Nx = Nx;
            single.Ny = Ny;
            single.Re = Re;
            single.nu = 1.0 / Re;
            single.U_lid = 0.0;
            single.U_inlet = U_inlet;
            single.t_end = tmax;
            single.dt = (equation_type == 9 || equation_type == 10) ? sys.getDouble("dt") : 0.01;
            single.is_steady = (equation_type == 8);
            single.fluid_mask = create_backward_step_mask(single.Nx, single.Ny, single.Lx, single.Ly);
        }
        tests.push_back(single);
    }

    // Запуск всех тестов по очереди
    for (const auto& tcase : tests) {
        run_single_test(tcase, equation_type, simp, piso, pimple, step_out);
    }

    std::cout << "\nAll tests finished.\n";
    return 0;
}