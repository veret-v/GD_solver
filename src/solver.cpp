#include "solver.h"

#include <vector>
#include <bits/stdc++.h>
#include <limits>

//RimanSolver1D::RimanSolver1D(Grid* grid) {

//}



double calc_sound_speed(const std::vector<double>& v_noncons, double g){
    return std::sqrt(g * v_noncons[P] / v_noncons[RHO]);
}

double calc_time_step(const std::vector<std::vector<double>>& v_cons, double dx, double cfl, double g, int fict_cells){
    double max_u = 0.0;
    const size_t Nx = v_cons.size();
    const size_t start = std::min(static_cast<size_t>(fict_cells), Nx);
    const size_t end = Nx > static_cast<size_t>(fict_cells) ? Nx - static_cast<size_t>(fict_cells) : 0;

    if (end <= start) {
        throw std::runtime_error("calc_time_step: not enough physical cells to advance solution");
    }

    for(size_t i = start; i < end; ++i){
        auto v_noncons = cons_to_noncons(v_cons[i], g);
        double c = calc_sound_speed(v_noncons, g);
        double u_val = v_cons[i][ru] / v_cons[i][r];
        double current_u = std::fabs(u_val) + c;
        if(current_u > max_u) {
            max_u = current_u;
        }
    }
    
    if (max_u <= std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("calc_time_step: maximum wave speed is zero or undefined");
    }
    
    return cfl * dx / max_u;
}

namespace {
constexpr double RHO_FLOOR = 1e-8;
constexpr double P_FLOOR = 1e-8;
constexpr double EPS = 1e-6;
constexpr int MAX_ITER_NUM = 20;
constexpr double P_MAX_RATIO = 2.0;
}

std::vector<double> cons_to_noncons(const std::vector<double>& v_cons, double g) {
    // �������� ������� ������
    
    if (v_cons.size() != M) {
        throw std::runtime_error("Invalid v_cons size in cons_to_noncons");
    }
    
    // �������� �� ������� ���������
    if (v_cons[r] <= 0.0) {
        throw std::runtime_error("Zero or negative density in cons_to_noncons");
    }
    
    // �������� ���������� ��������
    if (g <= 1.0) {
        throw std::runtime_error("Invalid gamma (g <= 1.0) in cons_to_noncons");
    }
    
    double rho = v_cons[r];
    double u = v_cons[ru] / rho;
    double p = (g - 1.0) * (v_cons[e] - 0.5 * v_cons[ru] * v_cons[ru] / rho);
    
    
    // �������� �� ������������� ��������
    if (p < 0.0) {
        std::cout << "WARNING: Negative pressure detected: " << p << std::endl;
        p = P_FLOOR; // ����� ������������� �������� ������ ��������������
    }
    
    return std::vector<double>{rho, u, p};
}

std::vector<double> noncons_to_cons(const std::vector<double>& v_noncons, double g) {
    // �������� ������� ������
    if (v_noncons.size() != M) {
        throw std::runtime_error("Invalid v_noncons size in noncons_to_cons");
    }
    
    // �������� ���������� ��������
    if (g <= 1.0) {
        throw std::runtime_error("Invalid gamma (g <= 1.0) in noncons_to_cons");
    }
    
    double rho = v_noncons[RHO];
    double u = v_noncons[U];
    double p = v_noncons[P];
    
    // �������� �� ������������� ���������
    if (rho < 0.0) {
        throw std::runtime_error("Negative density in noncons_to_cons");
    }
    
    double r_u = rho * u;
    double E = 0.5 * rho * u * u + p / (g - 1.0);
    
    return std::vector<double>{rho, r_u, E};
}

void enforce_physical_state(std::vector<double>& v_cons, double g) {
    if (v_cons.size() != M) {
        throw std::runtime_error("Invalid v_cons size in enforce_physical_state");
    }

    double rho = v_cons[r];
    double momentum = v_cons[ru];
    double energy = v_cons[e];

    if (rho < RHO_FLOOR) {
        double velocity = (rho > 0.0) ? momentum / rho : 0.0;
        rho = RHO_FLOOR;
        momentum = rho * velocity;
    }

    double velocity = momentum / rho;
    double kinetic = 0.5 * rho * velocity * velocity;
    double pressure = (g - 1.0) * (energy - kinetic);

    if (pressure < P_FLOOR) {
        pressure = P_FLOOR;
        energy = pressure / (g - 1.0) + kinetic;
    }

    v_cons[r] = rho;
    v_cons[ru] = rho * velocity;
    v_cons[e] = energy;
}


/* ��������� ����� */
std::vector<double> godunov_flux(const std::vector<double>& left_param, const std::vector<double>& right_param, double g){
    std::vector<double> v_noncons_l = cons_to_noncons(left_param, g);
    std::vector<double> v_noncons_r = cons_to_noncons(right_param, g);

    double r_l = v_noncons_l[RHO];
    double u_l = v_noncons_l[U];
    double p_l = v_noncons_l[P];

    double r_r = v_noncons_r[RHO];
    double u_r = v_noncons_r[U];
    double p_r = v_noncons_r[P];

    double c_l = calc_sound_speed(v_noncons_l, g);
    double c_r = calc_sound_speed(v_noncons_r, g);

    //double p_cont, v_cont;

    if ( 2. * (c_l + c_r) / (g - 1) <= u_r - u_l){
        std::cout << "error vacuum in godunov flux" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto [v_cont, p_cont] = calc_contact_pressure_velocity(v_noncons_l, v_noncons_r, c_l, c_r, g);

    
    std::vector<double> v_ncons = find_solution(v_noncons_l, v_noncons_r, v_cont, p_cont, c_l, c_r, 0.0, g);


    std::vector<double> flux = diff_flux_ncons(v_ncons, g);

    return flux;
    
}

std::vector<double> diff_flux_ncons(const std::vector<double>& v_ncons, double g) {
    double r_u = v_ncons[RHO] * v_ncons[U];
    double p_ru_2 = v_ncons[P] + (v_ncons[RHO] * v_ncons[U] * v_ncons[U]);
    double u_e_p = v_ncons[U] * (v_ncons[P] + 0.5 * v_ncons[RHO] * v_ncons[U] * v_ncons[U] + v_ncons[P] / (g - 1.0));
    std::vector<double> flux{r_u, p_ru_2, u_e_p};

    return flux;
}
 /* ��������� ������� */
std::vector<double> boundary(const std::vector<double>& v_cons, double boundary_type, double g) {
    // �������� ������� ������
    if (v_cons.size() != M) {
        throw std::runtime_error("Invalid v_cons size in boundary function");
    }
    
    
    std::vector<double> v_noncons = cons_to_noncons(v_cons, g);
    
    // �������� ���������� ��������������
    if (v_noncons.size() != M) {
        throw std::runtime_error("Invalid v_noncons size after conversion");
    }
    
    std::vector<double> bound(M, 0);
     // ����������� ������ ��� �������������
    if (boundary_type == 1) {
        // ���������� ������� (������) - �������� ������ ����
        bound[RHO] = v_noncons[RHO];
        bound[U] = -v_noncons[U];  // ��������� ��������
        bound[P] = v_noncons[P];
    }
    else if (boundary_type == 2) {
        // ��������� ����� (outflow) - ��� ��������� �����������
        bound[RHO] = v_noncons[RHO];
        bound[U] = v_noncons[U];  
        bound[P] = v_noncons[P];
    }
    else {
        // ��������� ������� �� ��������� ��� ����������� ���
        std::cout << "WARNING: Unknown boundary type " << boundary_type 
                  << ". Using free outflow." << std::endl;
        bound[RHO] = v_noncons[RHO];
        bound[U] = v_noncons[U];  
        bound[P] = v_noncons[P];
    }
    
    // ���������, ��� ��� ������ ���������
    if (bound.size() != M) {
        throw std::runtime_error("Boundary condition vector has wrong size");
    }
    std::vector<double> f = noncons_to_cons(bound, g);
    return f;
}



// ������� ��� ��������� ��������� ������� ����� ����
void set_sod_initial_conditions(std::vector<std::vector<double>>& u_prev, 
                               std::vector<std::vector<double>>& u_next,
                               int Nx, double x_min, double x_max, double rho_r, double u_r, double p_r, 
                                double rho_l, double u_l, double p_l, double g) {
    

    
    // ������ ������
    
    int fict_x = 1;
    int Nx_physical = Nx - 2 * fict_x;
    
    double dx = (x_max - x_min) / Nx_physical;
    
    //u_prev.resize(Nx, std::vector<double>(M, 0.0));
    //u_next.resize(Nx, std::vector<double>(M, 0.0));
   
    // ����� ������� (�������� �������)
    double x_mid = (x_min + x_max) / 2.0;
    
    for (int i = fict_x; i < Nx - fict_x; i++) {
        double x_cell = x_min + (i - fict_x + 0.5) * dx; // ���������� ������ ������
        
        if (x_cell < x_mid) {
            // ����� ����� - ������� �������� � ���������
            u_prev[i][r] = rho_l;  // ���������
            u_prev[i][ru] = rho_l * u_l;    // ��������
            u_prev[i][e] = 0.5 * rho_l * pow(u_l, 2.0) + p_l / (g - 1.0);    // ��������
        } else {
            // ������ ����� - ������ �������� � ���������
            u_prev[i][r] = rho_r; // ���������
            u_prev[i][ru] = rho_r * u_r;   // ��������
            u_prev[i][e] = 0.5 * rho_r * pow(u_r, 2.0) + p_r / (g - 1.0);    // ��������
        }
        
        // �������������� u_next ������ �� ����������
        u_next[i] = u_prev[i];
        enforce_physical_state(u_prev[i], g);
        enforce_physical_state(u_next[i], g);
        
    }
}

std::vector<double> find_solution(const std::vector<double>& v_noncons_l, const std::vector<double>& v_noncons_r,
                                  double v_cont, double p_cont, double c_l, double c_r, double S, double g){
    const double rl = v_noncons_l[RHO];
    const double vl = v_noncons_l[U];
    const double pl = v_noncons_l[P];

    const double rr = v_noncons_r[RHO];
    const double vr = v_noncons_r[U];
    const double pr = v_noncons_r[P];

    const double g1 = 0.5 * (g - 1.0) / g;
    const double g2 = 0.5 * (g + 1.0) / g;
    const double g3 = 2.0 * g / (g - 1.0);
    const double g4 = 2.0 / (g - 1.0);
    const double g5 = 2.0 / (g + 1.0);
    const double g6 = (g - 1.0) / (g + 1.0);
    const double g7 = 0.5 * (g - 1.0);

    double rho_star = 0.0;
    double vel_star = 0.0;
    double press_star = 0.0;

    if (S <= v_cont) {
        if (p_cont <= pl) {
            const double shl = vl - c_l;
            if (S <= shl) {
                rho_star = rl;
                vel_star = vl;
                press_star = pl;
            } else {
                const double p_ratio = p_cont / pl;
                const double cml = c_l * std::pow(p_ratio, g1);
                const double stl = v_cont - cml;
                if (S > stl) {
                    rho_star = rl * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (c_l + g7 * vl + S);
                    const double c = g5 * (c_l + g7 * (vl - S));
                    rho_star = rl * std::pow(c / c_l, g4);
                    press_star = pl * std::pow(c / c_l, g3);
                    vel_star = vel;
                }
            }
        } else {
            const double p_ratio = p_cont / pl;
            const double sl = vl - c_l * std::sqrt(g2 * p_ratio + g1);
            if (S <= sl) {
                rho_star = rl;
                vel_star = vl;
                press_star = pl;
            } else {
                rho_star = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        }
    } else {
        if (p_cont > pr) {
            const double p_ratio = p_cont / pr;
            const double sr = vr + c_r * std::sqrt(g2 * p_ratio + g1);
            if (S >= sr) {
                rho_star = rr;
                vel_star = vr;
                press_star = pr;
            } else {
                rho_star = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                vel_star = v_cont;
                press_star = p_cont;
            }
        } else {
            const double shr = vr + c_r;
            if (S >= shr) {
                rho_star = rr;
                vel_star = vr;
                press_star = pr;
            } else {
                const double p_ratio = p_cont / pr;
                const double cmr = c_r * std::pow(p_ratio, g1);
                const double str = v_cont + cmr;
                if (S <= str) {
                    rho_star = rr * std::pow(p_ratio, 1.0 / g);
                    vel_star = v_cont;
                    press_star = p_cont;
                } else {
                    const double vel = g5 * (-c_r + g7 * vr + S);
                    const double c = g5 * (c_r - g7 * (vr - S));
                    rho_star = rr * std::pow(c / c_r, g4);
                    press_star = pr * std::pow(c / c_r, g3);
                    vel_star = vel;
                }
            }
        }
    }

    std::vector<double> v_noncons(M);
    v_noncons[RHO] = rho_star;
    v_noncons[U] = vel_star;
    v_noncons[P] = press_star;
    return v_noncons;
}

double pressure_initial_guess(const std::vector<double>& v_ncons_l, const std::vector<double>& v_ncons_r, double cl, double cr, double g) {

    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
                         /* ���������� �������� */
    /* ��������� �����������, ������������ �� ����������� ������������ ��������������� �������
       � ����������� ���������� */
    double p_lin;
    double p_min, p_max;                /* ����������� � ������������ �������� ����� � ������ �� ������� */
    double p_ratio;                     /* ������� �� �������� ����� � ������ �� ������� */
    double p1, p2, g1, g2;              /* ��������������� ���������� ��� ������������� �������� */
    
    rl = v_ncons_l[RHO];
    vl = v_ncons_l[U];
    pl = v_ncons_l[P];
    rr = v_ncons_r[RHO];
    vr = v_ncons_r[U];
    pr = v_ncons_r[P];
    

    /* ��������� ����������� �� �������� ������
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = std::max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = std::min( pl, pr );
    p_max = std::max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( std::fabs( p_min - p_lin ) < EPS || std::fabs( p_max - p_lin ) < EPS ) ) ) {
        return p_lin;
    }

    if ( p_lin < p_min ) {
        g1 = 0.5 * ( g - 1.0 ) / g;
        return std::pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / std::pow( pl, g1 ) + cr / std::pow( pr, g1 ) ) ), 1.0 / g1 );
    }

    g1 = 2.0 / ( g + 1.0 );
    g2 = ( g - 1.0 ) / ( g + 1.0 );
    p1 = std::sqrt( g1 / rl / ( g2 * pl + p_lin ) );
    p2 = std::sqrt( g1 / rr / ( g2 * pr + p_lin ) );
    return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
    
}





// ��������������� ������� ��� �������� ���� ��������
std::tuple<double, double> calc_F_and_DF(double curr_press, const std::vector<double>& v_ncons, double c, double g) {
    double rho, v, p;
    double p_ratio, fg, q;
    double F, DF;

    rho = v_ncons[RHO];
    v = v_ncons[U];
    p = v_ncons[P];
    
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        // ����� ����������
        fg = 2.0 / (g - 1.0);
        F = fg * c * (std::pow(p_ratio, 1.0 / fg / g) - 1.0);
        DF = (1.0 / (rho * c)) * std::pow(p_ratio, -0.5 * (g + 1.0) / g);
    } else {
        // ������� �����
        q = std::sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
        F = (curr_press - p) / (c * rho * q);
        DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / (g * rho * c * std::pow(q, 3.0));
    }

    return std::make_tuple(F, DF);
}

// �������� �������, ������������ ������ (v_cont, p_cont)
std::tuple<double, double> calc_contact_pressure_velocity(
    const std::vector<double>& v_ncons_l, 
    const std::vector<double>& v_ncons_r,
    double cl, double cr, double g) 
{
    double vl, vr;
    double p_old;
    double fl, fr;
    double fld, frd;
    int iter_num = 0;
    double criteria;
    double p_cont, v_cont;

    vl = v_ncons_l[U];
    vr = v_ncons_r[U];
    
    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        printf("\ncalc_contact_pressure_velocity -> vacuum is generated\n");
        exit(EXIT_FAILURE);
    }

    // ������������, ��� pressure_initial_guess ���� ���������� ��� ����������
    p_old = pressure_initial_guess(v_ncons_l, v_ncons_r, cl, cr, g);
    if (p_old < 0.0) {
        printf("\ncalc_contact_pressure_velocity -> initial pressure guess is negative ");
        exit(EXIT_FAILURE);
    }
    
    // ������� ����������� ��������� ������� �������-�������
    do {
        // ���������� structured binding ��� ��������� ���� ��������
        auto [fl_temp, fld_temp] = calc_F_and_DF(p_old, v_ncons_l, cl, g);
        auto [fr_temp, frd_temp] = calc_F_and_DF(p_old, v_ncons_r, cr, g);
        
        fl = fl_temp;
        fld = fld_temp;
        fr = fr_temp;
        frd = frd_temp;
        
        p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
        criteria = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        iter_num++;
        
        if (iter_num > MAX_ITER_NUM) {
            printf("\ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value.\n");
            exit(EXIT_FAILURE);
        }
        if (p_cont < 0.0) {
            printf("\ncalc_contact_pressure_velocity -> pressure is negative.\n");
            exit(EXIT_FAILURE);            
        }
        p_old = p_cont;
    } while (criteria > EPS);

    // �������� ����������� �������
    v_cont = 0.5 * (vl + vr + fr - fl);

    return std::make_tuple(v_cont, p_cont);
}


