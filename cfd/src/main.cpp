#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "init.hpp"
#include "io.hpp"
#include "mesh_builder.hpp"
#include "schemes.hpp"

namespace cfd {

constexpr double kPi = 3.14159265358979323846;

struct Options {
    std::string mesh_type = "gmsh";
    std::string mesh_file = "meshes/naca0012_mesh.msh";
    std::string init = "freestream";
    std::string restart_file;
    std::string results_dir = "results";
    SpatialScheme scheme = SpatialScheme::Godunov;
    TimeScheme time_scheme = TimeScheme::RK2;

    double mach = 0.3;
    double alpha_deg = 0.0;
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
        << "  " << argv0 << " --mesh-type gmsh --mesh meshes/cylinder_channel.msh [options]\n\n"
        << "Options:\n"
        << "  --mesh <file>\n"
        << "  --mesh-type gmsh\n"
        << "  --mach <value>\n"
        << "  --alpha <deg>\n"
        << "  --gamma <value>\n"
        << "  --rho-inf <value>\n"
        << "  --p-inf <value>\n"
        << "  --cfl <value>\n"
        << "  --tend <value>\n"
        << "  --io <steps>\n"
        << "  --scheme godunov\n"
        << "  --time euler|rk2\n"
        << "  --init freestream|potential|restart\n"
        << "  --restart <file>\n"
        << "  --results <dir>\n"
        << "  --cyl-x <value> --cyl-y <value> --cyl-r <value>\n"
        << "  --ref-length <value>\n";
    std::exit(1);
}

double parse_double(const std::string &s)
{
    return std::stod(s);
}

int parse_int(const std::string &s)
{
    return std::stoi(s);
}

Options parse_args(int argc, char **argv)
{
    Options opt;
    opt.ref_length = 2.0 * opt.cyl_r;

    for(int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const std::string &name) -> std::string {
            if(i + 1 >= argc) throw std::runtime_error("Missing value after " + name);
            return argv[++i];
        };

        if(a == "--help" || a == "-h") {
            die_usage(argv[0]);
        } else if(a == "--mesh") {
            opt.mesh_file = next(a);
        } else if(a == "--mesh-type") {
            opt.mesh_type = next(a);
        } else if(a == "--mach") {
            opt.mach = parse_double(next(a));
        } else if(a == "--alpha") {
            opt.alpha_deg = parse_double(next(a));
        } else if(a == "--gamma") {
            opt.gamma = parse_double(next(a));
        } else if(a == "--rho-inf") {
            opt.rho_inf = parse_double(next(a));
        } else if(a == "--p-inf") {
            opt.p_inf = parse_double(next(a));
        } else if(a == "--cfl") {
            opt.cfl = parse_double(next(a));
        } else if(a == "--tend") {
            opt.tend = parse_double(next(a));
        } else if(a == "--io") {
            opt.io_each = parse_int(next(a));
        } else if(a == "--scheme") {
            const std::string s = to_lower(next(a));
            if(s != "godunov") throw std::runtime_error("Only --scheme godunov is implemented");
            opt.scheme = SpatialScheme::Godunov;
        } else if(a == "--time") {
            const std::string s = to_lower(next(a));
            if(s == "euler") opt.time_scheme = TimeScheme::Euler;
            else if(s == "rk2") opt.time_scheme = TimeScheme::RK2;
            else throw std::runtime_error("Unsupported time scheme: " + s);
        } else if(a == "--init") {
            opt.init = to_lower(next(a));
        } else if(a == "--restart") {
            opt.restart_file = next(a);
        } else if(a == "--results") {
            opt.results_dir = next(a);
        } else if(a == "--cyl-x") {
            opt.cyl_x = parse_double(next(a));
        } else if(a == "--cyl-y") {
            opt.cyl_y = parse_double(next(a));
        } else if(a == "--cyl-r") {
            opt.cyl_r = parse_double(next(a));
            opt.ref_length = 2.0 * opt.cyl_r;
        } else if(a == "--ref-length") {
            opt.ref_length = parse_double(next(a));
        } else {
            throw std::runtime_error("Unknown argument: " + a);
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

} // namespace cfd

int main(int argc, char **argv)
{
    using namespace cfd;

    try {
        const Options opt = parse_args(argc, argv);
        if(to_lower(opt.mesh_type) != "gmsh") {
            throw std::runtime_error("Only --mesh-type gmsh is implemented in this version");
        }

        std::filesystem::create_directories(opt.results_dir);

        FLOW.Mach = opt.mach;
        FLOW.alpha = opt.alpha_deg * kPi / 180.0;
        FLOW.gamma = opt.gamma;
        FLOW.rho_inf = opt.rho_inf;
        FLOW.p_inf = opt.p_inf;
        FLOW.init();

        Mesh mesh = load_gmsh(opt.mesh_file);
        initialize_field(mesh, opt);

        write_snap(mesh, opt.results_dir, 0, FLOW.gamma);

        std::ofstream aero(opt.results_dir + "/aero.dat");
        if(!aero) throw std::runtime_error("Failed to open aero.dat for writing");
        aero << "# t Cd Cl\n";

        double t = 0.0;
        int n = 0;
        int snap = 1;
        const double V_inf = std::sqrt(FLOW.u_inf * FLOW.u_inf + FLOW.v_inf * FLOW.v_inf);

        while(t < opt.tend - 1e-14) {
            double dt = compute_dt(mesh, opt.cfl, FLOW.gamma);
            dt = std::min(dt, opt.tend - t);

            step(mesh, dt, opt.scheme, opt.time_scheme, FLOW.gamma);
            t += dt;
            ++n;

            if(n % opt.io_each == 0 || t >= opt.tend - 1e-12) {
                write_snap(mesh, opt.results_dir, snap++, FLOW.gamma);
                const Aero ac = compute_aero(mesh, FLOW.rho_inf, V_inf, opt.ref_length, FLOW.gamma);
                aero << t << " " << ac.Cd << " " << ac.Cl << "\n";
                std::cout << "step=" << n << " t=" << t << " Cd=" << ac.Cd << " Cl=" << ac.Cl << "\n";
            }
        }

        write_vtk(mesh, opt.results_dir + "/out_final.vtk", FLOW.gamma);
        save_restart(mesh, opt.results_dir + "/restart.bin");
        write_cp(mesh, opt.results_dir + "/surface_cp.dat",
                 opt.cyl_x, opt.cyl_y, FLOW.rho_inf, V_inf, FLOW.p_inf, FLOW.gamma);

        std::cout << "Finished. Results written to " << opt.results_dir << "\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
