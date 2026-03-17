#ifndef MADER_H
#define MADER_H

#include <vector>

using Field3 = std::vector<std::vector<std::vector<double>>>;
using Mask2D = std::vector<std::vector<unsigned char>>;

struct MaderConfig {
    double visc = 0.5;
    bool step_enabled = false;
    double step_x_end = 0.0;
    double step_y_end = 0.0;
    double inflow_rho = 1.0;
    double inflow_u = 0.0;
    double inflow_v = 0.0;
    double inflow_p = 1.0;
};

Mask2D build_step_mask(double x_min, double y_min, double dx, double dy,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const MaderConfig& config);

void apply_step_mask(Field3& u,
                     const Mask2D& solid_mask,
                     const std::vector<double>& fill_cons);

double calc_time_step_masked(const Field3& v_cons,
                             const Mask2D& solid_mask,
                             double dx, double dy, double cfl, double g,
                             int fict_x, int fict_y);

void mader_step_2d(const Field3& u_prev,
                   Field3& u_next,
                   double dt,
                   double dx,
                   double dy,
                   double g,
                   int Nx,
                   int Ny,
                   int fict_x,
                   int fict_y,
                   const Mask2D& solid_mask,
                   const MaderConfig& config);

#endif
