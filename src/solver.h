#ifndef SOLVER
#define SOLVER

#include <vector>
#include "./grid.h"
#include "grid.h"



std::vector<double> godunov_flux(const std::vector<double>& left_param, const std::vector<double>& right_param, double g);

std::vector<double> cons_to_noncons(const std::vector<double>& v_cons, double g);

std::vector<double> noncons_to_cons(const std::vector<double>& v_noncons, double g);

double calc_time_step(const std::vector<std::vector<double>>& v_cons, double dx, double cfl, double g, int fict_cells);

double calc_sound_speed(const std::vector<double>& v_noncons, double g);

std::vector<double> boundary(const std::vector<double>& v_cons, double boundary, double g);

void enforce_physical_state(std::vector<double>& v_cons, double g);

void set_sod_initial_conditions(std::vector<std::vector<double>>& u_prev, 
                               std::vector<std::vector<double>>& u_next,
                               int Nx, double x_min, double x_max, double rho_r, double u_r, double p_r, 
                                double rho_l, double u_l, double p_l, double g); 

std::vector<double> find_solution(const std::vector<double>& v_noncons_l, const std::vector<double>& v_noncons_r,
                                  double v_cont, double p_cont, double c_l, double c_r, double S, double g);

std::tuple<double, double> calc_contact_pressure_velocity(
    const std::vector<double>& v_ncons_l, 
    const std::vector<double>& v_ncons_r,
    double cl, double cr, double g);
std::tuple<double, double> calc_F_and_DF(double curr_press, const std::vector<double>& v_ncons, double c, double g);  
double pressure_initial_guess(const std::vector<double>& v_ncons_l, const std::vector<double>& v_ncons_r, double cl, double cr, double g);
std::vector<double> diff_flux_ncons(const std::vector<double>& v_ncons, double g);

#endif