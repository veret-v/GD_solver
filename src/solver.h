#ifndef SOLVER
#define SOLVER

#define MAX_ITER_NUM 1000
#define EPS 1e-6
#define P_MAX_RATIO 2

#include "utils.h"
#include "grid.h"

#include <string>

class RimanSolver1D
{
private:
    Grid* grid;

    double g;
    double stop_time;
    double cl, cr; 
    double p_L, p_R, p_res;
    double rho_L, rho_R, rho_res;
    double u_L, u_R, u_res;
    double s;

    double p_cont;
    double v_cont;

    double fl, fr;
    double fld, frd;

    bool check_grid(const Grid &grid) ;

public:
    RimanSolver1D(
        Grid* grid, 
        double rho_L, double rho_R, 
        double p_L, double p_R, 
        double u_L, double u_R, 
        double g, double stop_time);

    void solve(const std::string &output_dir);
    void calc_sound_velocity();
    void calc_contact_pressure_velocity();
    void calc_F_and_DF(double curr_press);
    double pressure_initial_guess();
    void sample_solid_solution(Cell* cell);
};


#endif