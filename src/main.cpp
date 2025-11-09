#include <iostream>
#include <random>
#include <vector>
#include <filesystem>
#include <iomanip>
#include <fstream>

#include "./parser.h"
#include "./utils.h"
#include "./solver.h"
#include "./point.h"
#include "./grid.h"

int main(int argc, char** argv) {
    std::string system_ini = "../configs/initial.ini";
    std::string cases_ini = "../configs/SODA.ini";
    //std::string case_name = "case1";

    // ������� ���������������� ������
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
    std::string boundary_type_left = sys.getString("boundary_type_left");
    std::string boundary_type_right = sys.getString("boundary_type_right");
    std::string boundary_type_up = sys.getString("boundary_type_up");
    std::string boundary_type_down = sys.getString("boundary_type_down");
    std::string equation = sys.getString("equation");

    SectionedIniParser cases(cases_ini);
    double rho_L = cases.getDouble(case_name, "rho_L");
    double u_L = cases.getDouble(case_name, "u_L");
    double p_L = cases.getDouble(case_name, "p_L");
    double rho_R = cases.getDouble(case_name, "rho_R");
    double u_R = cases.getDouble(case_name, "u_R");
    double p_R = cases.getDouble(case_name, "p_R");

    double dx = (x_max - x_min) / Nx; 
    int Nx_with_fict_cells = Nx + 2 * fict_x;

    std::cout << "System ini: x_min=" << x_min << " x_max=" << x_max 
              << " y_min=" << y_min << " y_max=" << y_max 
              << " Nx=" << Nx << " Ny=" << Ny << " tmax=" << tmax
              << "\nCFL=" << cfl << " boundary_left=" << boundary_type_left 
              << " equation=" << equation << std::endl;
    std::cout << "Case [" << case_name << "]: "
              << "rho_L=" << rho_L << " u_L=" << u_L << " p_L=" << p_L
              << " rho_R=" << rho_R << " u_R=" << u_R << " p_R=" << p_R << std::endl;

    // ������������� ������
    std::vector<std::vector<double>> u_prev(Nx_with_fict_cells, std::vector<double>(M));
    std::vector<std::vector<double>> u_next(Nx_with_fict_cells, std::vector<double>(M));
    std::cout << "Initialization completed" << std::endl;
    std::cout << "u_prev size: " << u_prev.size() << std::endl;
    
    set_sod_initial_conditions(u_prev, u_next, Nx_with_fict_cells, x_min, x_max, 
                              rho_R, u_R, p_R, rho_L, u_L, p_L, g);

    

    // �������� ���� �� �������
    double curr_time = 0;
    int step = 0;
    
    while(curr_time < tmax) {
        // ��������� �������� ������
        for(int i = 0; i < fict_x; ++i) {
            u_prev[i] = boundary(u_prev[fict_x], 1, g);
            enforce_physical_state(u_prev[i], g);
            u_prev[Nx_with_fict_cells - 1 - i] = boundary(u_prev[Nx_with_fict_cells - fict_x - 1], 2, g);
            enforce_physical_state(u_prev[Nx_with_fict_cells - 1 - i], g);
        }

        double dt = calc_time_step(u_prev, dx, cfl, g, fict_x);

        std::vector<double> left_flux(M), right_flux(M);
        std::vector<double> bound(M);

        // ���������� ������� � ���������� �������
        for(int i = fict_x; i < Nx_with_fict_cells - fict_x; ++i) {
            if(i == fict_x) {
                bound = boundary(u_prev[i], 1, g);
                left_flux = godunov_flux(bound, u_prev[i], g);
            } else {
                left_flux = godunov_flux(u_prev[i-1], u_prev[i], g);
            }

            if(i == Nx_with_fict_cells - fict_x - 1) {
                bound = boundary(u_prev[i], 2, g);
                right_flux = godunov_flux(u_prev[i], bound, g);
            } else {
                right_flux = godunov_flux(u_prev[i], u_prev[i+1], g);
            }

            for(int j = 0; j < M; ++j) {
                u_next[i][j] = u_prev[i][j] - (dt / dx) * (right_flux[j] - left_flux[j]);
            }
            enforce_physical_state(u_next[i], g);
        }

        // ����������� ������ ������ ��������� �� ���������� ���������
        for(int i = 0; i < fict_x; ++i) {
            u_next[i] = boundary(u_next[fict_x], 1, g);
            enforce_physical_state(u_next[i], g);
            u_next[Nx_with_fict_cells - 1 - i] = boundary(u_next[Nx_with_fict_cells - fict_x - 1], 2, g);
            enforce_physical_state(u_next[Nx_with_fict_cells - 1 - i], g);
        }

        std::swap(u_prev, u_next);
        curr_time += dt;
        step++;

        // ����� ��������� ������ 100 �����
        if(step % 100 == 0) {
            std::cout << "Step: " << step << ", Time: " << curr_time 
                      << ", dt: " << dt << std::endl;
        }
    }

    // ������ ����������� � ����
    std::string output_dir = "../output/";
    std::filesystem::create_directories(output_dir);
    
    std::ofstream fout(output_dir + "final_results.csv");
    fout << "x,rho,u,p\n";
    
    for(int i = fict_x; i < Nx_with_fict_cells - fict_x; i++) {
        double x = x_min + (i - fict_x + 0.5) * dx;
        std::vector<double> prim = cons_to_noncons(u_prev[i], g);
        fout << x << "," << prim[RHO] << "," << prim[U] << "," << prim[P] << "\n";
    }
    fout.close();

    std::cout << "Calculation completed!" << std::endl;
    std::cout << "Final time: " << curr_time << std::endl;
    std::cout << "Total steps: " << step << std::endl;
    std::cout << "Results saved to: " << output_dir + "final_results.csv" << std::endl;

    return 0;
}


