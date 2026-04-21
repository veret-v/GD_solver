#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "init.hpp"
#include "io.hpp"
#include "mesh_builder.hpp"
#include "schemes.hpp"
#include <mpi.h>
#include "mpi_parallel.hpp"
namespace cfd {

constexpr double kPi = 3.14159265358979323846;

struct Options {
    std::string case_type = "cylinder"; // "cylinder" или "airfoil"
    std::string mesh_type = "gmsh";
    std::string mesh_file = "meshes/cylinder_channel.msh";
    std::string init = "freestream";
    std::string restart_file;
    std::string results_dir = "results";
    SpatialScheme scheme = SpatialScheme::Godunov;
    TimeScheme time_scheme = TimeScheme::RK2;

    double mach = 0.3;
    double alpha_deg = 5.0;
    double gamma = 1.4;
    double rho_inf = 1.0;
    double p_inf = 1.0 / 1.4;
    double cfl = 0.4;
    double tend = 5.0;
    int io_each = 50;

    double cyl_x = 1.5;
    double cyl_y = 1.0;
    double cyl_r = 0.25;
    double ref_length = 0.5;
};

[[noreturn]] void die_usage(const char *argv0)
{
    std::cerr
        << "Usage:\n"
        << "  " << argv0 << " --case cylinder|airfoil --mesh meshes/grid.msh [options]\n\n"
        << "Options:\n"
        << "  --case cylinder|airfoil\n"
        << "  --mesh <file>\n"
        << "  --mach <value>\n"
        << "  --alpha <deg>\n"
        << "  --cfl <value>\n"
        << "  --tend <value>\n"
        << "  --io <steps>\n"
        << "  --scheme godunov|kolgan\n"
        << "  --time euler|rk2\n"
        << "  --init freestream|potential|restart\n"
        << "  --results <dir>\n"
        << "  --cyl-x <value> --cyl-y <value> --cyl-r <value>\n"
        << "  --ref-length <value>\n";
    std::exit(1);
}

double parse_double(const std::string &s) { return std::stod(s); }
int parse_int(const std::string &s) { return std::stoi(s); }

Options parse_args(int argc, char **argv)
{
    Options opt;
    opt.ref_length = 2.0 * opt.cyl_r; // Default to cylinder diameter

    for(int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const std::string &name) -> std::string {
            if(i + 1 >= argc) throw std::runtime_error("Missing value after " + name);
            return argv[++i];
        };

        if(a == "--help" || a == "-h") die_usage(argv[0]);
        else if(a == "--case") opt.case_type = to_lower(next(a));
        else if(a == "--mesh") opt.mesh_file = next(a);
        else if(a == "--mesh-type") opt.mesh_type = next(a);
        else if(a == "--mach") opt.mach = parse_double(next(a));
        else if(a == "--alpha") opt.alpha_deg = parse_double(next(a));
        else if(a == "--gamma") opt.gamma = parse_double(next(a));
        else if(a == "--rho-inf") opt.rho_inf = parse_double(next(a));
        else if(a == "--p-inf") opt.p_inf = parse_double(next(a));
        else if(a == "--cfl") opt.cfl = parse_double(next(a));
        else if(a == "--tend") opt.tend = parse_double(next(a));
        else if(a == "--io") opt.io_each = parse_int(next(a));
        else if(a == "--scheme") {
            const std::string s = to_lower(next(a));
            if(s == "godunov") opt.scheme = SpatialScheme::Godunov;
            else if(s == "kolgan") opt.scheme = SpatialScheme::Kolgan;
            else throw std::runtime_error("Unsupported scheme: " + s);
        } else if(a == "--time") {
            const std::string s = to_lower(next(a));
            if(s == "euler") opt.time_scheme = TimeScheme::Euler;
            else if(s == "rk2") opt.time_scheme = TimeScheme::RK2;
            else throw std::runtime_error("Unsupported time scheme: " + s);
        } else if(a == "--init") opt.init = to_lower(next(a));
        else if(a == "--restart") opt.restart_file = next(a);
        else if(a == "--results") opt.results_dir = next(a);
        else if(a == "--cyl-x") opt.cyl_x = parse_double(next(a));
        else if(a == "--cyl-y") opt.cyl_y = parse_double(next(a));
        else if(a == "--cyl-r") {
            opt.cyl_r = parse_double(next(a));
            opt.ref_length = 2.0 * opt.cyl_r;
        } else if(a == "--ref-length") opt.ref_length = parse_double(next(a));
        else throw std::runtime_error("Unknown argument: " + a);
    }

    if (opt.case_type == "airfoil") {
        if (opt.init == "potential") {
            std::cerr << "Warning: Potential flow is for cylinders. Falling back to freestream init.\n";
            opt.init = "freestream";
        }
    }

    return opt;
}

void initialize_field(Mesh &mesh, const Options &opt)
{
    if(opt.init == "freestream") {
        init_freestream(mesh);
    } else if(opt.init == "potential") {
        init_potential(mesh, opt.cyl_x, opt.cyl_y, opt.cyl_r, FLOW.gamma);
    } else if(opt.init == "restart") {
        if(opt.restart_file.empty()) throw std::runtime_error("--init restart requires --restart");
        load_restart(mesh, opt.restart_file);
    } else {
        throw std::runtime_error("Unsupported init mode: " + opt.init);
    }
}

// Запись Cp для профиля (вдоль оси x)
void write_airfoil_cp(const Mesh &m, const std::string &filename, int rank, double rho_inf, double V_inf, double p_inf, double g)
{
    std::string ranked_name = filename + "_r" + std::to_string(rank) + ".dat";
    std::ofstream f(ranked_name);
    const double q = 0.5 * rho_inf * V_inf * V_inf;
    std::vector<std::pair<double, double>> data;
    
    for(const Face &f : m.faces) {
        if(f.bc != Face::BC::Wall) continue;
        const Primitive W = conservative_to_primitive(m.cells[f.left].U, g);
        const double cp = (W.p - p_inf) / q;
        data.emplace_back(f.mx, cp); // x-координата вместо угла
    }

    std::sort(data.begin(), data.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    
    if(!f) throw std::runtime_error("Failed to open Cp file: " + filename);
    f << "# x Cp\n";
    for(const auto &[x, cp] : data) f << x << " " << cp << "\n";
}

} // namespace cfd

int main(int argc, char **argv)
{
    using namespace cfd;
    MPI_Init(&argc, &argv);

    try {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        const Options opt = parse_args(argc, argv);
        if (rank == 0) std::filesystem::create_directories(opt.results_dir);

        FLOW.Mach = opt.mach;
        FLOW.alpha = opt.alpha_deg * kPi / 180.0;
        FLOW.gamma = opt.gamma;
        FLOW.rho_inf = opt.rho_inf;
        FLOW.p_inf = opt.p_inf;
        FLOW.init();

        // 1. Каждый читает сетку
        Mesh global_mesh = load_gmsh(opt.mesh_file);
        std::vector<int> part(global_mesh.nc(), 0);

        // 2. Нулевой ранг делает декомпозицию RCB
        if (rank == 0) {
            std::vector<int> ids(global_mesh.nc());
            for(int i = 0; i < global_mesh.nc(); ++i) ids[i] = i;
            partition_rcb(ids, size, 0, global_mesh.cells, part);
            std::cout << "RCB Partitioning completed.\n";
        }

        // 3. Рассылаем массив принадлежности
        MPI_Bcast(part.data(), part.size(), MPI_INT, 0, MPI_COMM_WORLD);

        // 4. Каждый строит локальную часть и настраивает Halo
        Halo halo;
        Mesh mesh = build_local_mesh(global_mesh, part, rank, halo);
        
        initialize_field(mesh, opt);
        write_snap(mesh, opt.results_dir, 0, FLOW.gamma);

        std::ofstream aero;
        if (rank == 0) {
            aero.open(opt.results_dir + "/aero.dat");
            aero << "# t Cd Cl\n";
        }

        double t = 0.0;
        int n = 0, snap = 1;
        const double V_inf = std::sqrt(FLOW.u_inf * FLOW.u_inf + FLOW.v_inf * FLOW.v_inf);
        double start_time = MPI_Wtime();
        while(t < opt.tend - 1e-14) {
            double dt = compute_dt(mesh, opt.cfl, FLOW.gamma);
            dt = std::min(dt, opt.tend - t);

            // Передаем halo
            step(mesh, halo, dt, opt.scheme, opt.time_scheme, FLOW.gamma);
            t += dt;
            ++n;

            if(n % opt.io_each == 0 || t >= opt.tend - 1e-12) {
                write_snap(mesh, opt.results_dir, snap++, FLOW.gamma);
                
                const Aero ac_raw = compute_aero(mesh, FLOW.rho_inf, V_inf, opt.ref_length, FLOW.gamma);
                
                if (rank == 0) {
                    double Cd =  ac_raw.Cd * std::cos(FLOW.alpha) + ac_raw.Cl * std::sin(FLOW.alpha);
                    double Cl = -ac_raw.Cd * std::sin(FLOW.alpha) + ac_raw.Cl * std::cos(FLOW.alpha);
                    aero << t << " " << Cd << " " << Cl << "\n";
                    std::cout << "step=" << n << " t=" << t << " Cd=" << Cd << " Cl=" << Cl << "\n";
                }
            }
        }

        if (rank == 0) std::cout << "Finished. Results written to " << opt.results_dir << "\n";
        
        // Очистка типов
        MPI_Type_free(&halo.cell_mpi_type);

        // После завершения цикла и записи финальных файлов:
        double end_time = MPI_Wtime();
        double duration = end_time - start_time;

        // Выводим результат только от лица главного процесса (rank 0)
        if (rank == 0) {
            std::cout << "\n========================================" << std::endl;
            std::cout << "End!" << std::endl;
            std::cout << "Time: " << std::fixed << std::setprecision(3) 
                      << duration << " s." << std::endl;
            
            // Если итераций было много, можно посчитать скорость:
            // std::cout << "Среднее время на шаг: " << duration / n << " сек." << std::endl;
            std::cout << "========================================\n" << std::endl;
        }

    } catch(const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Finalize();
    return 0;
}
