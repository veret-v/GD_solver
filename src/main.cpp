// main.cpp — последовательная реализация SIMPLE / PISO / PIMPLE
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
void apply_boundary_conditions(std::vector<std::vector<double>>& u,
                               std::vector<std::vector<double>>& v,
                               int Nx, int Ny, double U_lid)
{
    // u имеет размер (Nx+1) x Ny
    // Левая стенка (i = 0) – неподвижна
    for (int j = 0; j < Ny; ++j) u[0][j] = 0.0;
    // Правая стенка (i = Nx) – неподвижна
    for (int j = 0; j < Ny; ++j) u[Nx][j] = 0.0;
    // Верхняя крышка (j = Ny-1) – движется со скоростью U_lid
    for (int i = 0; i <= Nx; ++i) u[i][Ny-1] = U_lid;
    // Нижняя стенка для u не требует явной установки,
    // так как граничное условие уже учтено в коэффициентах.

    // v имеет размер Nx x (Ny+1)
    // Нижняя стенка (j = 0) – неподвижна
    for (int i = 0; i < Nx; ++i) v[i][0] = 0.0;
    // Верхняя стенка (j = Ny) – неподвижна
    for (int i = 0; i < Nx; ++i) v[i][Ny] = 0.0;
    // Боковые стенки для v (i = 0 и i = Nx-1) – неподвижны
    for (int j = 1; j < Ny; ++j) {
        v[0][j] = 0.0;
        v[Nx-1][j] = 0.0;
    }
}

int main(int argc, char** argv) {
    // Чтение конфигурации
    std::string system_ini = "../configs/input.ini";
    std::string cases_ini = "../configs/SODA.ini";  // или другой файл с начальными условиями

    IniParser sys(system_ini);
    std::string case_name = sys.getString("case_name");
    double x_min = sys.getDouble("x_min");
    double x_max = sys.getDouble("x_max");
    double y_min = sys.getDouble("y_min");
    double y_max = sys.getDouble("y_max");
    int Nx = sys.getInt("Nx");
    int Ny = sys.getInt("Ny");
    double tmax = sys.getDouble("tmax");
    double cfl = sys.getDouble("cfl");   // для несжимаемых методов не используется напрямую
    double Re = sys.getDouble("Re");      // число Рейнольдса
    double U_lid = sys.getDouble("U_lid"); // скорость крышки (для каверны)
    int equation_type = sys.getInt("equation_type"); // 8=SIMPLE, 9=PISO, 10=PIMPLE

    // Параметры метода
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

    // Параметры вывода
    int step_out = sys.getInt("step_out");
    std::string output_dir = "../output/steps/";
    std::filesystem::create_directories(output_dir);

    // Физические параметры
    double rho = 1.0;
    double nu = 1.0 / Re;
    double Lx = x_max - x_min;
    double Ly = y_max - y_min;
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    // Инициализация staggered массивов
    std::vector<std::vector<double>> u(Nx+1, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> v(Nx, std::vector<double>(Ny+1, 0.0));
    std::vector<std::vector<double>> p(Nx, std::vector<double>(Ny, 0.0));

    // Начальные условия (для каверны: везде 0, кроме верхней границы u = U_lid)
    for (int i = 0; i <= Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            u[i][j] = 0.0;
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j <= Ny; ++j)
            v[i][j] = 0.0;
    // Верхняя крышка
    apply_boundary_conditions(u, v, Nx, Ny, U_lid);

    // Функция для сохранения текущего состояния в CSV
    auto save_state = [&](int step, double time, const std::string& suffix = "") {
        std::ostringstream filename;
        filename << output_dir << "step_" << std::setw(5) << std::setfill('0') << step
                 << "_t_" << std::fixed << std::setprecision(4) << time << suffix << ".csv";
        std::ofstream fout(filename.str());
        fout << "x,y,u,v,p\n";
        // Вычисляем скорость в центрах ячеек для вывода
        for (int i = 0; i < Nx; ++i) {
            double x = x_min + (i + 0.5) * dx;
            for (int j = 0; j < Ny; ++j) {
                double y = y_min + (j + 0.5) * dy;
                double u_c = 0.5 * (u[i][j] + u[i+1][j]);
                double v_c = 0.5 * (v[i][j] + v[i][j+1]);
                fout << x << "," << y << ","
                     << u_c << "," << v_c << "," << p[i][j] << "\n";
            }
        }
        fout.close();
    };

    // Сохраняем начальное состояние
    save_state(0, 0.0, "_initial");

    auto start_time = std::chrono::high_resolution_clock::now();

    // =============== SIMPLE ===============
    if (equation_type == 8) {
        std::cout << "Starting SIMPLE, Nx=" << Nx << ", Ny=" << Ny << ", Re=" << Re << "\n";
        for (int iter = 1; iter <= simp.max_iter; ++iter) {
            // Сохраняем старые значения для проверки сходимости
            auto u_old = u;
            auto v_old = v;

            // 1. Коэффициенты и предиктор u*
            std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
            compute_momentum_coefficients_u(u, v, p, dx, dy, nu, rho, 0.0, true,
                                            aP_u, aE_u, aW_u, aN_u, aS_u, H_u, U_lid);
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);
            std::vector<std::vector<double>> u_star = u;
            for (int i = 1; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    double H = H_u[i][j];
                    double dp = (p[i][j] - p[i-1][j]) * dy;
                    u_star[i][j] = (1.0 - simp.alpha_u) * u[i][j] +
                                    (simp.alpha_u / aP_u[i][j]) * (H - dp);
                }
            }

            // 2. Коэффициенты и предиктор v*
            std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
            compute_momentum_coefficients_v(u, v, p, dx, dy, nu, rho, 0.0, true,
                                            aP_v, aE_v, aW_v, aN_v, aS_v, H_v);
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);
            std::vector<std::vector<double>> v_star = v;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 1; j < Ny; ++j) {
                    double H = H_v[i][j];
                    double dp = (p[i][j] - p[i][j-1]) * dx;
                    v_star[i][j] = (1.0 - simp.alpha_u) * v[i][j] +
                                    (simp.alpha_u / aP_v[i][j]) * (H - dp);
                }
            }

            // 3. Уравнение Пуассона для p'
            std::vector<std::vector<double>> p_corr(Nx, std::vector<double>(Ny, 0.0));
            solve_pressure_correction(u_star, v_star, aP_u, aP_v, dx, dy, p_corr);

            // 4. Коррекция
            correct_velocity_pressure(u, v, p, u_star, v_star, p_corr, aP_u, aP_v,
                                      dx, dy, simp.alpha_p);
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);

            // Проверка сходимости
            double max_du = 0.0;
            for (int i = 1; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    max_du = std::max(max_du, std::abs(u[i][j] - u_old[i][j]));
            for (int i = 0; i < Nx; ++i)
                for (int j = 1; j < Ny; ++j)
                    max_du = std::max(max_du, std::abs(v[i][j] - v_old[i][j]));

            if (iter % 100 == 0 || max_du < simp.tol) {
                std::cout << "SIMPLE iter " << iter << ", max_du = " << max_du << "\n";
                if (iter % step_out == 0 || max_du < simp.tol)
                    save_state(iter, iter, "");
            }
            if (max_du < simp.tol) {
                std::cout << "SIMPLE converged at iteration " << iter << "\n";
                break;
            }
        }
    }

    // =============== PISO ===============
    else if (equation_type == 9) {
        std::cout << "Starting PISO, Nx=" << Nx << ", Ny=" << Ny << ", Re=" << Re
                  << ", dt=" << piso.dt << ", t_end=" << piso.t_end << "\n";
        double t = 0.0;
        int step = 0;
        save_state(0, 0.0, "");
        while (t < piso.t_end) {
            auto u_old = u;
            auto v_old = v;
            auto p_old = p;

            // Предиктор (без релаксации)
            std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
            compute_momentum_coefficients_u(u_old, v_old, p_old, dx, dy, nu, rho,
                                            piso.dt, false, aP_u, aE_u, aW_u, aN_u, aS_u, H_u, U_lid);
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);
            std::vector<std::vector<double>> u_star = u_old;
            for (int i = 1; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    double H = H_u[i][j];
                    double dp = (p_old[i][j] - p_old[i-1][j]) * dy;
                    u_star[i][j] = (H - dp) / aP_u[i][j];
                }
            }

            std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
            compute_momentum_coefficients_v(u_old, v_old, p_old, dx, dy, nu, rho,
                                            piso.dt, false, aP_v, aE_v, aW_v, aN_v, aS_v, H_v);
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);
            std::vector<std::vector<double>> v_star = v_old;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 1; j < Ny; ++j) {
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
                std::vector<std::vector<double>> p_corr(Nx, std::vector<double>(Ny, 0.0));
                solve_pressure_correction(u, v, aP_u, aP_v, dx, dy, p_corr);
                correct_velocity_pressure(u, v, p, u, v, p_corr, aP_u, aP_v, dx, dy, 1.0);
            }
            apply_boundary_conditions(u, v, Nx, Ny, U_lid);

            t += piso.dt;
            step++;

            if (step % step_out == 0) {
                std::cout << "PISO step " << step << ", t = " << t << "\n";
                save_state(step, t, "");
            }
        }
    }

    // =============== PIMPLE ===============
    else if (equation_type == 10) {
        std::cout << "Starting PIMPLE, Nx=" << Nx << ", Ny=" << Ny << ", Re=" << Re
                  << ", dt=" << pimple.dt << ", t_end=" << pimple.t_end
                  << ", nOuter=" << pimple.n_outer << ", nCorr=" << pimple.n_correctors << "\n";
        double t = 0.0;
        int step = 0;
        save_state(0, 0.0, "");
        while (t < pimple.t_end) {            
            auto u_old = u;
            auto v_old = v;
            auto p_old = p;

            for (int outer = 0; outer < pimple.n_outer; ++outer) {
                double alpha_u_cur = (outer == pimple.n_outer - 1) ? 1.0 : pimple.alpha_u;
                double alpha_p_cur = (outer == pimple.n_outer - 1) ? 1.0 : pimple.alpha_p;

                // Коэффициенты u
                std::vector<std::vector<double>> aP_u, aE_u, aW_u, aN_u, aS_u, H_u;
                compute_momentum_coefficients_u(u, v, p, dx, dy, nu, rho,
                                                pimple.dt, false, aP_u, aE_u, aW_u, aN_u, aS_u, H_u, U_lid);
                apply_boundary_conditions(u, v, Nx, Ny, U_lid);
                std::vector<std::vector<double>> u_star = u;
                for (int i = 1; i < Nx; ++i) {
                    for (int j = 0; j < Ny; ++j) {
                        double H = H_u[i][j];
                        double dp = (p[i][j] - p[i-1][j]) * dy;
                        u_star[i][j] = (1.0 - alpha_u_cur) * u[i][j] +
                                        (alpha_u_cur / aP_u[i][j]) * (H - dp);
                        // Сохраняем эффективную диагональ
                        aP_u[i][j] /= alpha_u_cur;
                    }
                }

                // Коэффициенты v
                std::vector<std::vector<double>> aP_v, aE_v, aW_v, aN_v, aS_v, H_v;
                compute_momentum_coefficients_v(u, v, p, dx, dy, nu, rho,
                                                pimple.dt, false, aP_v, aE_v, aW_v, aN_v, aS_v, H_v);
                apply_boundary_conditions(u, v, Nx, Ny, U_lid);
                std::vector<std::vector<double>> v_star = v;
                for (int i = 0; i < Nx; ++i) {
                    for (int j = 1; j < Ny; ++j) {
                        double H = H_v[i][j];
                        double dp = (p[i][j] - p[i][j-1]) * dx;
                        v_star[i][j] = (1.0 - alpha_u_cur) * v[i][j] +
                                        (alpha_u_cur / aP_v[i][j]) * (H - dp);
                        aP_v[i][j] /= alpha_u_cur;
                    }
                }

                // Внутренние коррекции
                std::vector<std::vector<double>> u_corr = u_star;
                std::vector<std::vector<double>> v_corr = v_star;
                std::vector<std::vector<double>> p_corr_field = p;
                for (int corr = 0; corr < pimple.n_correctors; ++corr) {
                    std::vector<std::vector<double>> p_corr(Nx, std::vector<double>(Ny, 0.0));
                    solve_pressure_correction(u_corr, v_corr, aP_u, aP_v, dx, dy, p_corr);
                    correct_velocity_pressure(u_corr, v_corr, p_corr_field,
                                              u_corr, v_corr, p_corr,
                                              aP_u, aP_v, dx, dy, alpha_p_cur);
                }
                apply_boundary_conditions(u, v, Nx, Ny, U_lid);
                u = u_corr;
                v = v_corr;
                p = p_corr_field;
            }

           

            t += pimple.dt;
            step++;

            if (step % step_out == 0) {
                std::cout << "PIMPLE step " << step << ", t = " << t << "\n";
                save_state(step, t, "");
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Calculation completed in " << duration.count() / 1000.0 << " seconds.\n";

    // Сохраняем финальное состояние
    save_state(9999, (equation_type == 8) ? simp.max_iter : tmax, "_final");

    return 0;
}
/*
#include <iostream>
#include <random>
#include <vector>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include "./parser.h"
#include "./utils.h"
#include "./solver.h"
#include "./point.h"
#include "./grid.h"
#include <script.h>
#include <mpi.h>
int main(int argc, char** argv) {
// Изменено на единый файл input.ini
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string system_ini = "../configs/input.ini";
    std::string cases_ini = "../configs/SODA.ini";

    // чтение системных параметров из объединенного файла
    IniParser sys(system_ini);
    std::string case_name = sys.getString("case_name");
    double x_min = sys.getDouble("x_min");
    double x_max = sys.getDouble("x_max");
    double y_min = sys.getDouble("y_min");
    double y_max = sys.getDouble("y_max");
    int Nx = sys.getInt("Nx");
    int Ny = sys.getInt("Ny");
    int fict_x = sys.getInt("fict_x");
    int fict_y = sys.getInt("fict_y");
    double tmax = sys.getDouble("tmax");
    double cfl = sys.getDouble("cfl");
    double g = sys.getDouble("g");
    double dt_out = sys.getDouble("dt_out");
    int step_num = sys.getInt("step_num");
    int px = sys.getInt("px");
    int py = sys.getInt("py");
    std::string boundary_type_left = sys.getString("boundary_type_left");
    std::string boundary_type_right = sys.getString("boundary_type_right");
    std::string boundary_type_up = sys.getString("boundary_type_up");
    std::string boundary_type_down = sys.getString("boundary_type_down");
    int equation_type = sys.getInt("equation_type");

    SectionedIniParser cases(cases_ini);
    double rho_L = cases.getDouble(case_name, "rho_L");
    double u_L = cases.getDouble(case_name, "u_L");
    double v_L = cases.getDouble(case_name, "v_L");
    double p_L = cases.getDouble(case_name, "p_L");
    double rho_R = cases.getDouble(case_name, "rho_R");
    double u_R = cases.getDouble(case_name, "u_R");
    double v_R = cases.getDouble(case_name, "v_R");
    double p_R = cases.getDouble(case_name, "p_R");

    int dims[2] = {px, py};
    MPI_Dims_create(size, 2, dims); // dims[0] - по X, dims[1] - по Y

    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int left_rank, right_rank, down_rank, up_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &left_rank, &right_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &down_rank, &up_rank);

    // 2. Вычисление РАБОЧИХ границ для текущего процесса в глобальных координатах
    int i_start_phys, i_end_phys;
    int j_start_phys, j_end_phys;

    get_subdomain_bounds(Nx, dims[0], coords[0], i_start_phys, i_end_phys);
    get_subdomain_bounds(Ny, dims[1], coords[1], j_start_phys, j_end_phys);
    
    // Сдвигаем индексы на размер фиктивных ячеек
    int i_start = i_start_phys + fict_x;
    int i_end   = i_end_phys + fict_x;
    int j_start = j_start_phys + fict_y;
    int j_end   = j_end_phys + fict_y;
    
    double dx = (x_max - x_min) / Nx; 
    double dy = (y_max - y_min) / Ny;
    int Nx_with_fict_cells = Nx + 2 * fict_x;
    int Ny_with_fict_cells = Ny + 2 * fict_y;

 

    // инициализация сетки
    std::vector<std::vector<std::vector<double>>> u_prev(Nx_with_fict_cells,  std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> u_next(Nx_with_fict_cells,  std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    parseAndInitialize(u_prev, "C:\\solvver\\GD_solver\\configs\\input.ini", fict_x, fict_y, g);
    parseAndInitialize(u_next, "C:\\solvver\\GD_solver\\configs\\input.ini", fict_x, fict_y, g);
    //set_sod_initial_conditions(u_prev, u_next, Nx, x_min, x_max, Ny, y_min, y_max,
      //                        rho_R, u_R, v_R, p_R, rho_L, u_L, v_L, p_L, g, fict_x, fict_y);

    auto parse_boundary_code = [](const std::string& value) -> int {
        try {
            int code = std::stoi(value);
            return (code <= 0) ? 2 : code;
        } catch (...) {
            return 2; // default to outflow
        }
    };

    const int left_bc_code = parse_boundary_code(boundary_type_left);
    const int right_bc_code = parse_boundary_code(boundary_type_right);
    const int up_bc_code = parse_boundary_code(boundary_type_up);
    const int down_bc_code = parse_boundary_code(boundary_type_down);

    std::string steps_dir = "../output/steps/";
 if (rank == 0) {
        std::filesystem::create_directories(steps_dir);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Сборка и сохранение начальных условий
    gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    if (rank == 0) {
        std::ofstream fout_initial(steps_dir + "step_0_initial.csv");
        fout_initial << "x,y,rho,u,v,p\n";
        for(int i = fict_x; i < Nx_with_fict_cells - fict_x; i++) {
            double x = x_min + (i - fict_x + 0.5) * dx;
            for(int j = fict_y; j < Ny_with_fict_cells - fict_y; j++) {
                double y = y_min + (j - fict_y + 0.5) * dy;
                std::vector<double> prim = cons_to_noncons(u_prev[i][j], g);
                fout_initial << x << "," << y << "," << prim[RHO] << "," << prim[U] << "," << prim[V] << "," << prim[P] << "\n";
            }
        }
        fout_initial.close();
    }

    // основной цикл по времени
    double curr_time = 0;
    int step = 0;

// Массивы для хранения реконструкции
    // X-direction structures
    std::vector<std::vector<std::vector<double>>> left_face_x(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> right_face_x(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> delta_x(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    
    // Y-direction structures (NEW)
    // left_face_y соответствует "нижней" грани (j-1/2), right_face_y - "верхней" (j+1/2) в локальных координатах
    std::vector<std::vector<std::vector<double>>> left_face_y(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> right_face_y(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> delta_y(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    
    // Массивы реконструкции для Родионова (Corrector step arrays)
    std::vector<std::vector<std::vector<double>>> rec_x_L(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_x_R(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_y_L(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> rec_y_R(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    std::vector<std::vector<std::vector<double>>> u_half(Nx_with_fict_cells, std::vector<std::vector<double>>(Ny_with_fict_cells, std::vector<double>(M)));
    auto start = std::chrono::high_resolution_clock::now();
    while(curr_time < tmax) {


        exchange_halos_global(u_prev, cart_comm, fict_x, fict_y, i_start, i_end, j_start, j_end, left_rank, right_rank, down_rank, up_rank);
     
        // 1. Установка граничных условий на текущем слое
        // 1. Граничные условия (Ghost Cells)
        // Лево/Право
        for(int j = 0; j < Ny_with_fict_cells; ++j) {
            for(int k = 0; k < fict_x; ++k) {
                u_prev[k][j] = boundary(u_prev[fict_x][j], left_bc_code, g, 0); // Axis 0 (X)
                u_prev[Nx_with_fict_cells - 1 - k][j] = boundary(u_prev[Nx_with_fict_cells - fict_x - 1][j], right_bc_code , g, 0);
            }
        }
        // Верх/Низ
        for(int i = 0; i < Nx_with_fict_cells; ++i) {
            for(int k = 0; k < fict_y; ++k) {
                u_prev[i][k] = boundary(u_prev[i][fict_y], down_bc_code , g, 1); // Axis 1 (Y)
                u_prev[i][Ny_with_fict_cells - 1 - k] = boundary(u_prev[i][Ny_with_fict_cells - fict_y - 1], up_bc_code , g, 1);
            }
        }

        double dt = calc_time_step(u_prev, dx, dy, cfl, g, fict_x);
        if ((curr_time + dt) > tmax){
            dt = tmax - curr_time; // Исправлена логика обрезки шага
        }
       

        // 2. Физические граничные условия (применяются только если нет соседа)
        if (left_rank == MPI_PROC_NULL) {
            for(int j = std::max(0, j_start - 1); j <= std::min(Ny_with_fict_cells - 1, j_end); ++j)
                for(int k = 0; k < fict_x; ++k)
                    u_prev[k][j] = boundary(u_prev[fict_x][j], left_bc_code, g, 0); 
        }
        if (right_rank == MPI_PROC_NULL) {
            for(int j = std::max(0, j_start - 1); j <= std::min(Ny_with_fict_cells - 1, j_end); ++j)
                for(int k = 0; k < fict_x; ++k)
                    u_prev[Nx_with_fict_cells - 1 - k][j] = boundary(u_prev[Nx_with_fict_cells - fict_x - 1][j], right_bc_code , g, 0);
        }
        if (up_rank == MPI_PROC_NULL) {
            for(int i = std::max(0, i_start - 1); i <= std::min(Nx_with_fict_cells - 1, i_end); ++i)
                for(int k = 0; k < fict_y; ++k)
                    u_prev[i][Ny_with_fict_cells - 1 - k] = boundary(u_prev[i][Ny_with_fict_cells - fict_y - 1], up_bc_code, g, 1);
        }
        if (down_rank == MPI_PROC_NULL) {
            for(int i = std::max(0, i_start - 1); i <= std::min(Nx_with_fict_cells - 1, i_end); ++i)
                for(int k = 0; k < fict_y; ++k)
                    u_prev[i][k] = boundary(u_prev[i][fict_y], down_bc_code, g, 1);
        }
        // 3. Синхронизация шага по времени
        double local_dt = 1e10; // Вычисляем ТОЛЬКО по своему поддомену
        for(int i = i_start; i < i_end; ++i) {
            for(int j = j_start; j < j_end; ++j) {
                auto prim = cons_to_noncons(u_prev[i][j], g);
                double c = calc_sound_speed(prim, g);
                double sig_x = std::abs(prim[U]) + c;
                double sig_y = std::abs(prim[V]) + c;
                double inv_dt = (sig_x / dx) + (sig_y / dy);
                local_dt = std::min(local_dt, cfl / inv_dt);
            }
        }

        double dt;
        MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
        
        // --- ПРЕДВАРИТЕЛЬНЫЕ РАСЧЕТЫ (Вне цикла по ячейкам) ---
        // Расширенные границы для реконструкции (чтобы получить грани для расчета потоков)
        int rec_i_start = std::max(0, i_start - 1);
        int rec_i_end   = std::min(Nx_with_fict_cells - 1, i_end);
        int rec_j_start = std::max(0, j_start - 1);
        int rec_j_end   = std::min(Ny_with_fict_cells - 1, j_end);

        // TYPE 1: KOLGAN (2D Extension)
        if(equation_type == 1) {
            // 1. Reconstruct along X (for each row j)
            
            for(int j = rec_j_start; j <= rec_j_end; ++j) {
                for(int i = rec_i_start; i <= rec_i_end  ; ++i) {
                     if(i == 0) {
                        Kolgan(u_prev[i][j], u_prev[i][j], u_prev[i+1][j], left_face_x[i][j], right_face_x[i][j], g);
                    } 
                    else if(i == Nx_with_fict_cells - 1) {
                         Kolgan(u_prev[i-1][j], u_prev[i][j], 
                                u_prev[i][j], left_face_x[i][j], 
                                right_face_x[i][j], g);
                     }
                     else{
                        Kolgan(u_prev[i-1][j], u_prev[i][j], u_prev[i+1][j], 
                                left_face_x[i][j], right_face_x[i][j], g);
                     }
                }
            }
            // 2. Reconstruct along Y (for each column i)
            for(int i = rec_i_start; i <= rec_i_end; ++i) {
                for(int j = rec_j_start; j <= rec_j_end; ++j) {
                    if (j == 0){
                        Kolgan(u_prev[i][j], u_prev[i][j], u_prev[i][j+1], left_face_y[i][j], right_face_y[i][j], g);
                        }
                    else if (j == Ny_with_fict_cells - 1){
                        Kolgan(u_prev[i][j-1], 
                                u_prev[i][j], 
                                u_prev[i][j], left_face_y[i][j], right_face_y[i][j], g);
                        }
                    else{
                        Kolgan(u_prev[i][j-1], u_prev[i][j], u_prev[i][j+1], left_face_y[i][j], right_face_y[i][j], g);
                        }
                }
            }
        }

        // TYPE 2: RODIONOV (2D Extension)
        if(equation_type == 2) {
            // 1. Расчет дельт (наклонов) по X
            for(int j = rec_j_start; j <= rec_j_end; ++j) {
                for(int i = rec_i_start; i <= rec_i_end; ++i) {
                    if(i == 0) {
                        Kolgan_for_Rodionov(u_prev[0][j], u_prev[0][j], u_prev[1][j], delta_x[i][j], g);
                    } 
                    else if( i == Nx_with_fict_cells - 1){
                        Kolgan_for_Rodionov(u_prev[Nx_with_fict_cells - 2][j], 
                                                u_prev[Nx_with_fict_cells - 1][j], u_prev[Nx_with_fict_cells - 1][j],
                                                delta_x[i][j], g);
                    }
                    else{
                        Kolgan_for_Rodionov(u_prev[i-1][j], u_prev[i][j], u_prev[i+1][j], delta_x[i][j], g);
                    }
                }
            }
            // 2. Расчет дельт (наклонов) по Y
            for(int i = rec_i_start; i <= rec_i_end; ++i) {
                for(int j = rec_j_start; j <= rec_j_end; ++j) {
                    if(j == 0){
                        Kolgan_for_Rodionov(u_prev[i][0], u_prev[i][0], u_prev[i][1], delta_y[i][j], g);
                    }
                    else if(j == Ny_with_fict_cells - 1){
                        Kolgan_for_Rodionov(u_prev[i][Ny_with_fict_cells - 2],
                                            u_prev[i][Ny_with_fict_cells - 1], 
                                            u_prev[i][Ny_with_fict_cells - 1], delta_y[i][j], g);
                    }
                    else{
                        Kolgan_for_Rodionov(u_prev[i][j-1], u_prev[i][j], u_prev[i][j+1], delta_y[i][j], g);
                    }
                }
            }

            // 3. Шаг Предиктор (расчет u_half) - теперь полноценный 2D
            Rodionov(u_prev, u_half, delta_x, delta_y, dt, dx, dy, g, 
                     Nx_with_fict_cells, Ny_with_fict_cells, fict_x, fict_y,
                     left_bc_code, right_bc_code, up_bc_code, down_bc_code);

            // 4. Шаг Корректор: Реконструкция на гранях
            for (int i = rec_i_start; i <= rec_i_end; i++) {
                for (int j = rec_j_start; j <= rec_j_end; j++) {
                    for (int var = 0; var < M; var++) {
                        double u_avg = 0.5 * (u_prev[i][j][var] + u_half[i][j][var]);
                        
                        // По X: rec_x_L это "левое" значение внутри ячейки (для потока слева это right state)
                        rec_x_L[i][j][var] = u_avg - 0.5 * delta_x[i][j][var];
                        rec_x_R[i][j][var] = u_avg + 0.5 * delta_x[i][j][var];
                        
                        // По Y:
                        rec_y_L[i][j][var] = u_avg - 0.5 * delta_y[i][j][var];
                        rec_y_R[i][j][var] = u_avg + 0.5 * delta_y[i][j][var];
                    }
                    enforce_physical_state(rec_x_L[i][j], g);
                    enforce_physical_state(rec_x_R[i][j], g);
                    enforce_physical_state(rec_y_L[i][j], g);
                    enforce_physical_state(rec_y_R[i][j], g);
                }
            }
        }
           

        // --- ОСНОВНОЙ ЦИКЛ ПО ЯЧЕЙКАМ (Расчет потоков и обновление) ---
        for(int i = i_start; i < i_end; ++i) {
            for(int j = j_start; j < j_end; ++j){
            std::vector<double> left_flux(M), right_flux(M);
            std::vector<double> up_flux(M), down_flux(M);

            // Выбор схемы для расчета потоков
            if(equation_type == 0){ // Godunov
                left_flux = godunov_flux_x(u_prev[i-1][j], u_prev[i][j], g);
                right_flux = godunov_flux_x(u_prev[i][j], u_prev[i+1][j], g);
                up_flux = godunov_flux_y(u_prev[i][j], u_prev[i][j+1], g);
                down_flux = godunov_flux_y(u_prev[i][j-1], u_prev[i][j], g);
            }
            else if(equation_type == 1){ // Kolgan
                // Используем заранее рассчитанные грани 
                left_flux = godunov_flux_x(right_face_x[i-1][j], left_face_x[i][j], g);
                right_flux = godunov_flux_x(right_face_x[i][j], left_face_x[i+1][j], g);
                up_flux = godunov_flux_y(right_face_y[i][j], left_face_y[i][j+1], g);
                down_flux = godunov_flux_y(right_face_y[i][j-1], left_face_y[i][j], g);
            }
            else if(equation_type == 2){ // Rodionov
                // Используем заранее рассчитанные left/right 
                left_flux = godunov_flux_x(rec_x_R[i-1][j], rec_x_L[i][j], g);
                right_flux = godunov_flux_x(rec_x_R[i][j], rec_x_L[i+1][j], g);
                up_flux = godunov_flux_y(rec_y_R[i][j], rec_y_L[i][j+1], g);
                down_flux = godunov_flux_y(rec_y_R[i][j-1], rec_y_L[i][j], g);
            }
            else if (equation_type == 3){ // HLL
                left_flux = hll_flux_new(u_prev[i-1][j], u_prev[i][j], g, 0);
                right_flux = hll_flux_new(u_prev[i][j], u_prev[i+1][j], g, 0);
                up_flux = hll_flux_new(u_prev[i][j], u_prev[i][j+1], g, 1);
                down_flux = hll_flux_new(u_prev[i][j-1], u_prev[i][j], g, 1);
            }
            else if(equation_type == 4){ // HLLC
                left_flux = hllc_flux_new(u_prev[i-1][j], u_prev[i][j], g, 0);
                right_flux = hllc_flux_new(u_prev[i][j], u_prev[i+1][j], g, 0);
                up_flux = hllc_flux_new(u_prev[i][j], u_prev[i][j+1], g, 1);
                down_flux = hllc_flux_new(u_prev[i][j-1], u_prev[i][j], g, 1);
            }
            else if(equation_type == 5){ // Rusanov 
                left_flux = rusanov_2d(u_prev[i-1][j], u_prev[i][j], g, 0);
                right_flux = rusanov_2d(u_prev[i][j], u_prev[i+1][j], g, 0);
                up_flux = rusanov_2d(u_prev[i][j], u_prev[i][j+1], g, 1);
                down_flux = rusanov_2d(u_prev[i][j-1], u_prev[i][j], g, 1);
            }
            else if(equation_type == 6){ // Osher
                left_flux = osher_flux_2d(u_prev[i-1][j], u_prev[i][j], g, 0);
                right_flux = osher_flux_2d(u_prev[i][j], u_prev[i+1][j], g, 0);
                up_flux = osher_flux_2d(u_prev[i][j], u_prev[i][j+1], g, 1);
                down_flux = osher_flux_2d(u_prev[i][j-1], u_prev[i][j], g, 1);
            }
            else if(equation_type == 7){ // Roe
                left_flux = roe_flux_2d(u_prev[i-1][j], u_prev[i][j], g, 0);
                right_flux = roe_flux_2d(u_prev[i][j], u_prev[i+1][j], g, 0);
                up_flux = roe_flux_2d(u_prev[i][j], u_prev[i][j+1], g, 1);
                down_flux = roe_flux_2d(u_prev[i][j-1], u_prev[i][j], g, 1);
            }

            // Обновление решения u_next
            for(int k = 0; k < M; ++k) {
                // Стандартный конечно-объемный апдейт
                u_next[i][j][k] = u_prev[i][j][k] - (dt / dx) * (right_flux[k] - left_flux[k])
                                    - (dt / dy) * (up_flux[k] - down_flux[k]);
            }
            
            enforce_physical_state(u_next[i][j], g);
            }
        }
          

        std::swap(u_prev, u_next);
        curr_time += dt;
        step++;
    
    
    
    
    
    


// Вывод и сохранение с использованием MPI-сборки
        if(step % step_num == 0) {
            if (rank == 0) {
                std::cout << "Step: " << step << ", Time: " << curr_time << ", dt: " << dt << std::endl;
            }
            
            gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);

            if (rank == 0) {
                std::ostringstream filename;
                filename << steps_dir << "step_" << step << "_time_" << std::fixed << std::setprecision(6) << curr_time << ".csv";
                
                std::vector<double> left_prim{rho_L, u_L, v_L, p_L};
                std::vector<double> right_prim{rho_R, u_R, v_R, p_R};
                auto analytic_solution = compute_analytic_solution_2d(x_min, x_max, y_min, y_max, Nx, Ny, fict_x, fict_y,
                                                                 curr_time, left_prim, right_prim, g);
                
                std::ofstream fout_step(filename.str());
                fout_step << "x,y,rho,u,v,p,rho_exact,u_exact,v_exact,p_exact\n";
                for(int i = fict_x; i < Nx_with_fict_cells - fict_x; i++) {
                    double x = x_min + (i - fict_x + 0.5) * dx;
                    for(int j = fict_y; j < Ny_with_fict_cells - fict_y; j++) {
                        double y = y_min + (j - fict_y + 0.5) * dy;
                        std::vector<double> prim = cons_to_noncons(u_prev[i][j], g);
                        std::vector<double> prim_exact = analytic_solution[i][j];
                    
                        fout_step << x << "," << y << "," << prim[RHO] << "," << prim[U] << "," << prim[V] << "," <<prim[P] << ","
                             << prim_exact[0] << "," << prim_exact[1] << "," << prim_exact[2] << "," << prim_exact[3] <<"\n";
                    }
                }
                fout_step.close();
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start);
    
    // Сборка и сохранение финальных результатов
    gather_to_root(u_prev, Nx, Ny, fict_x, fict_y, i_start, i_end, j_start, j_end, rank, size, cart_comm, dims);
    
    if (rank == 0) {
        std::string output_dir = "../output/steps/";
        std::filesystem::create_directories(output_dir);
        
        std::vector<double> left_prim{rho_L, u_L, v_L, p_L};
        std::vector<double> right_prim{rho_R, u_R, v_R, p_R};
        auto analytic_final = compute_analytic_solution_2d(
        x_min, x_max, y_min, y_max, Nx, Ny, fict_x, fict_y, curr_time, left_prim, right_prim, g);
        
        std::ofstream fout(output_dir + "final_results.csv");
        fout << "x,y,rho,u,v,p,rho_exact,u_exact,v_exact,p_exact\n";
        for(int i = fict_x; i < Nx_with_fict_cells - fict_x; i++) {
            double x = x_min + (i - fict_x + 0.5) * dx;
            for(int j = fict_y; j < Ny_with_fict_cells - fict_y; j++){
                double y = y_min + (j - fict_y + 0.5) * dy;
                std::vector<double> prim = cons_to_noncons(u_prev[i][j], g);
                std::vector<double> prim_exact = analytic_final[i][j];
            
                fout << x << "," << y << "," << prim[RHO] << "," << prim[U] << "," << prim[V] << "," << prim[P] << ","
                     << prim_exact[0] << "," << prim_exact[1] << "," << prim_exact[2] << "," << prim_exact[3] <<"\n";
            }
        }
        fout.close();

        double seconds = duration.count() / 1000.0;
        std::cout << "==========================================" << std::endl;
        std::cout << "Calculation completed! Final time: " << curr_time << std::endl;
        std::cout << "Total execution time: " << std::fixed << std::setprecision(3) 
                  << seconds << " seconds (" << duration.count() << " milliseconds)" << std::endl;
        std::cout << "==========================================" << std::endl;
    }

    MPI_Finalize();
    return 0;

}
*/