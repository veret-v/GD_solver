#include <iostream>
#include <random>
#include <vector>
#include <filesystem>
#include <iomanip>

#include "./parser.h"
#include "./utils.h"
#include "./solver.h"
#include "./point.h"
#include "./grid.h"

int main(int argc, char** argv) {
    std::string system_ini = "../configs/initial.ini";
    std::string cases_ini = "../configs/SODA.ini";

    std::string case_name  = "case1";
    // if (argc > 1) system_ini = argv[1];
    // if (argc > 2) cases_ini  = argv[2];
    // if (argc > 3) case_name  = argv[3];

    IniParser sys(system_ini);
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
    std::string boundary_type_left  = sys.getString("boundary_type_left");
    std::string boundary_type_right = sys.getString("boundary_type_right");
    std::string boundary_type_up    = sys.getString("boundary_type_up");
    std::string boundary_type_down  = sys.getString("boundary_type_down");
    std::string equation = sys.getString("equation");

    SectionedIniParser cases(cases_ini);
    double rho_L = cases.getDouble(case_name, "rho_L");
    double u_L   = cases.getDouble(case_name, "u_L");
    double p_L   = cases.getDouble(case_name, "p_L");
    double rho_R = cases.getDouble(case_name, "rho_R");
    double u_R   = cases.getDouble(case_name, "u_R");
    double p_R   = cases.getDouble(case_name, "p_R");

    std::cout << "System ini: x_min=" << x_min << " x_max=" << x_max << " y_min=" << y_min << " y_max=" << y_max << " Nx=" << Nx <<  "Ny=" << Ny << " tmax=" << tmax
              << "\nCFL=" << cfl << " boundary=" << boundary_type_left << " equation=" << equation << std::endl;
    std::cout << "Case [" << case_name << "]: "
              << "rho_L=" << rho_L << " u_L=" << u_L << " p_L=" << p_L
              << " rho_R=" << rho_R << " u_R=" << u_R << " p_R=" << p_R << std::endl;


    std::string output_dir = "../output/data.csv";
    Grid task(fict_x, fict_y, Nx, Ny, x_min, x_max, y_min, y_max, boundary_type_left, boundary_type_right, boundary_type_down, boundary_type_up);
    task.set_values(rho_L, rho_R, p_L, p_R, Point(u_L, 0), Point(u_R, 0));
    task.WriteCSV(output_dir);
}
