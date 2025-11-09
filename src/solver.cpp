#include "solver.h"


RimanSolver1D::RimanSolver1D(
    Grid* grid, 
    double rho_L, double rho_R, 
    double p_L, double p_R, 
    double u_L, double u_R, 
    double g, double stop_time
){
    if (check_grid(*grid)) this -> grid = grid;
    this -> rho_L = rho_L;
    this -> rho_R = rho_R;
    this -> p_L = p_L;
    this -> p_R = p_R;
    this -> u_L = u_L;
    this -> u_R = u_R;
    this -> g = g;
    this -> stop_time = stop_time;
}

bool RimanSolver1D::check_grid(const Grid &grid) 
{
    if(grid.cells_num_y > 1)
        return 0;
    return 1;
}

void RimanSolver1D::calc_sound_velocity() 
{
    this -> cl = sqrt( g * p_L / rho_L);
    this -> cr = sqrt( g * p_R / rho_R );

}

void RimanSolver1D::calc_contact_pressure_velocity() {
    double p_old;       
    double fl, fr;      
    double fld, frd;    
    int iter_num = 0;   
    double criteria;   
    double g;           

    
    if ( 2.0 * ( cl + cr ) / ( g - 1.0 ) <= u_R - u_L ) {
        printf( "\ncalc_contact_pressure_velocity -> vacuum is generated\n" );
        exit( EXIT_FAILURE );
    }

    p_old = pressure_initial_guess();
    if ( p_old < 0.0 ) {
        printf( "\ncalc_contact_pressure_velocity -> initial pressure guess is negative " );
        exit( EXIT_FAILURE );
    }
    
    do {
        calc_F_and_DF(p_old);
        p_cont = p_old - ( fl + fr + u_R - u_L ) / ( fld + frd );
        criteria = 2.0 * fabs( (p_cont - p_old ) / (p_cont + p_old ) );
        iter_num++;
        if ( iter_num > MAX_ITER_NUM ) {
            printf( "\ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value.\n" );
            exit( EXIT_FAILURE );
        }
        if ( p_cont < 0.0 ) {
            printf( "\ncalc_contact_pressure_velocity -> pressure is negative.\n" );
            exit( EXIT_FAILURE );            
        }
        p_old = p_cont;
    } while ( criteria > EPS );

    v_cont = 0.5 * ( u_L + u_R + fr - fl );
}

void RimanSolver1D::calc_F_and_DF(double curr_press)
{
    double p_ratio, fg, q;         

    p_ratio = curr_press / p_L;
    if ( curr_press <= p_L ) {
        fg = 2.0 / ( g - 1.0 );
        fl = fg * cl * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        fld = ( 1.0 / rho_L / cl ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        fl = ( curr_press - p_L ) / cl / rho_L / q;
        fld = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / rho_L / cl / pow( q, 3.0 );
    }

    p_ratio = curr_press / p_R;
    if ( curr_press <= p_R ) {
        fg = 2.0 / ( g - 1.0 );
        fr = fg * cr * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        frd = ( 1.0 / rho_R / cr ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        fr = ( curr_press - p_R ) / cr / rho_R / q;
        frd = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / rho_R / cr / pow( q, 3.0 );
    }
}


double RimanSolver1D::pressure_initial_guess() 
{
    double p_lin;
    double p_min, p_max;               
    double p_ratio;                     
    double p1, p2, g1, g2;              
    
    p_lin = max(0.0, 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L) * (rho_L + rho_R) * ( cl + cr ) );
    p_min = min(p_L, p_R);
    p_max = max(p_L, p_R);
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < EPS || fabs( p_max - p_lin ) < EPS ) ) ) {
        return p_lin;
    } else {
        if ( p_lin < p_min ) {
            g1 = 0.5 * ( g - 1.0 ) / g;
            return pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * (u_R - u_L) ) / ( cl / pow(p_L, g1 ) + cr / pow(p_R, g1 ) ) ), 1.0 / g1 );
        } else {
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = sqrt( g1 / rho_L / ( g2 * p_L + p_lin ) );
            p2 = sqrt( g1 / rho_R / ( g2 * p_R + p_lin ) );
            return ( p1 * p_L + p2 * p_R - ( u_R - u_L ) ) / ( p1 + p2 );
        }
    }
}

void RimanSolver1D::sample_solid_solution(Cell* cell) 
{              
    double g1, g2, g3, g4, g5, g6, g7;      

    double shl, stl;        
    double sl, sr;              

    double shr, str;        

    double cml, cmr;        
    double c;               
    double p_ratio;
    double r, v, p;         
    
    g1 = 0.5 * ( g - 1.0 ) / g;
    g2 = 0.5 * ( g + 1.0 ) / g;
    g3 = 2.0 * g / ( g - 1.0 );
    g4 = 2.0 / ( g - 1.0 );
    g5 = 2.0 / ( g + 1.0 );
    g6 = ( g - 1.0 ) / ( g + 1.0 );
    g7 = 0.5 * ( g - 1.0 );

    if ( s <= v_cont ) {
        if ( p_cont <= p_L ) {
            shl = u_L - cl;
            if ( s <= shl ) {
                r = rho_L;
                v = u_L;
                p = p_L;
            }
            else {
                cml = cl * pow( p_cont / p_L, g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    r = rho_L * pow( p_cont / p_L, 1.0 / g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    v = g5 * ( cl + g7 * u_L + s );
                    c = g5 * ( cl + g7 * ( u_L - s ) );
                    r = rho_L * pow( c / cl, g4 );
                    p = p_L * pow( c / cl, g3 );
                }
            }
        }
        else {
            p_ratio = p_cont / p_L;
            sl = u_L - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                r = rho_L;
                v = u_L;
                p = p_L;
            }
            else {
                r = rho_L * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        if ( p_cont > p_R ) {
            p_ratio = p_cont / p_R;
            sr = u_R + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                r = rho_R;
                v = u_R;
                p = p_R;
            }
            else {
                r = rho_R * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            shr = u_R + cr;
            if ( s >= shr ) {
                r = rho_R;
                v = u_R;
                p = p_R;
            }
            else {
               cmr = cr * pow( p_cont / p_R, g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   r = rho_R * pow( p_cont / p_R, 1.0 / g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    v = g5 * ( - cr + g7 * u_R + s );
                    c = g5 * ( cr - g7 * ( u_R - s ) );
                    r = rho_R * pow( c / cr, g4 );
                    p = p_R * pow( c / cr, g3 );
               }
            }
        }
    }
     
    cell->rho = r;
    cell->u.x = v;
    cell->p = p;
}

void RimanSolver1D::solve(const std::string &output_dir)
{
    for (auto cell : grid -> grid) {
        s = cell.center.x / stop_time;
        sample_solid_solution(&cell);
    }
    (*grid).WriteCSV(output_dir);
}