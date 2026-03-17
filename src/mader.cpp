#include "mader.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "grid.h"
#include "solver.h"

namespace {

using Field2 = std::vector<std::vector<double>>;

constexpr double EPS = 1e-8;

struct FaceState {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double E = 0.0;
    bool blocked = true;
    bool physical = false;
    int i = -1;
    int j = -1;
};

bool is_physical(int i, int j, int Nx, int Ny, int fict_x, int fict_y) {
    return i >= fict_x && i < Nx - fict_x && j >= fict_y && j < Ny - fict_y;
}

bool is_fluid_cell(int i, int j, int Nx, int Ny, int fict_x, int fict_y, const Mask2D& solid_mask) {
    return is_physical(i, j, Nx, Ny, fict_x, fict_y) && solid_mask[i][j] == 0;
}

std::vector<double> make_inflow_cons(const MaderConfig& config, double g) {
    return noncons_to_cons(std::vector<double>{
        config.inflow_rho,
        config.inflow_u,
        config.inflow_v,
        config.inflow_p
    }, g);
}

Field2 make_field2(int Nx, int Ny, double value = 0.0) {
    return Field2(Nx, std::vector<double>(Ny, value));
}

std::vector<unsigned char> build_step_rows(const Mask2D& solid_mask,
                                           int Nx, int Ny,
                                           int fict_x, int fict_y) {
    std::vector<unsigned char> step_rows(Ny, 0);
    for (int j = fict_y; j < Ny - fict_y; ++j) {
        for (int i = fict_x; i < Nx - fict_x; ++i) {
            if (solid_mask[i][j]) {
                step_rows[j] = 1;
                break;
            }
        }
    }
    return step_rows;
}

void apply_mader_boundaries(Field3& u,
                            int Nx, int Ny, int fict_x, int fict_y,
                            const Mask2D& solid_mask,
                            const MaderConfig& config,
                            double g) {
    const std::vector<double> inflow_cons = make_inflow_cons(config, g);
    const std::vector<unsigned char> step_rows = build_step_rows(solid_mask, Nx, Ny, fict_x, fict_y);

    for (int j = 0; j < Ny; ++j) {
        const bool wall_left = (j >= fict_y && j < Ny - fict_y && step_rows[j] != 0);
        for (int k = 0; k < fict_x; ++k) {
            if (wall_left) {
                u[k][j] = boundary(u[fict_x][j], 1, g, 0);
            } else {
                u[k][j] = inflow_cons;
            }
            u[Nx - 1 - k][j] = boundary(u[Nx - fict_x - 1][j], 2, g, 0);
        }
    }

    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < fict_y; ++k) {
            u[i][k] = boundary(u[i][fict_y], 1, g, 1);
            u[i][Ny - 1 - k] = boundary(u[i][Ny - fict_y - 1], 1, g, 1);
        }
    }

    apply_step_mask(u, solid_mask, inflow_cons);
}

std::vector<double> sample_cons(const Field3& u,
                                int i, int j,
                                int ni, int nj,
                                int Nx, int Ny,
                                int fict_x, int fict_y,
                                const Mask2D& solid_mask,
                                const std::vector<unsigned char>& step_rows,
                                const std::vector<double>& inflow_cons,
                                double g) {
    if (ni < fict_x) {
        if (step_rows[j]) {
            return boundary(u[i][j], 1, g, 0);
        }
        return inflow_cons;
    }
    if (ni >= Nx - fict_x) {
        return boundary(u[i][j], 2, g, 0);
    }
    if (nj < fict_y || nj >= Ny - fict_y) {
        return boundary(u[i][j], 1, g, 1);
    }
    if (solid_mask[ni][nj]) {
        const int axis = (ni != i) ? 0 : 1;
        return boundary(u[i][j], 1, g, axis);
    }
    return u[ni][nj];
}

double face_u_sum_left(int i, int j,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const Mask2D& solid_mask,
                       const std::vector<unsigned char>& step_rows,
                       const Field2& u_old,
                       const Field2& u_tilde,
                       const MaderConfig& config) {
    if (i == fict_x) {
        return step_rows[j] ? 0.0 : 2.0 * config.inflow_u;
    }
    if (!is_fluid_cell(i - 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return 0.0;
    }
    return u_old[i - 1][j] + u_tilde[i - 1][j];
}

double face_u_sum_right(int i, int j,
                        int Nx, int Ny, int fict_x, int fict_y,
                        const Mask2D& solid_mask,
                        const Field2& u_old,
                        const Field2& u_tilde) {
    if (i == Nx - fict_x - 1) {
        return u_old[i][j] + u_tilde[i][j];
    }
    if (!is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return 0.0;
    }
    return u_old[i + 1][j] + u_tilde[i + 1][j];
}

double face_v_sum_down(int i, int j,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const Mask2D& solid_mask,
                       const Field2& v_old,
                       const Field2& v_tilde) {
    if (j == fict_y) {
        return 0.0;
    }
    if (!is_fluid_cell(i, j - 1, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return 0.0;
    }
    return v_old[i][j - 1] + v_tilde[i][j - 1];
}

double face_v_sum_up(int i, int j,
                     int Nx, int Ny, int fict_x, int fict_y,
                     const Mask2D& solid_mask,
                     const Field2& v_old,
                     const Field2& v_tilde) {
    if (j == Ny - fict_y - 1) {
        return 0.0;
    }
    if (!is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return 0.0;
    }
    return v_old[i][j + 1] + v_tilde[i][j + 1];
}

double clamp_fraction(double value) {
    return std::max(-0.999, std::min(0.999, value));
}

void accumulate_transport(Field2& DM,
                          Field2& DE,
                          Field2& DPU,
                          Field2& DPV,
                          const FaceState& donor,
                          const FaceState& acceptor,
                          double amount) {
    if (amount <= 0.0) {
        return;
    }

    if (donor.physical) {
        DM[donor.i][donor.j] -= amount;
        DE[donor.i][donor.j] -= donor.E * amount;
        DPU[donor.i][donor.j] -= donor.u * amount;
        DPV[donor.i][donor.j] -= donor.v * amount;
    }
    if (acceptor.physical) {
        DM[acceptor.i][acceptor.j] += amount;
        DE[acceptor.i][acceptor.j] += donor.E * amount;
        DPU[acceptor.i][acceptor.j] += donor.u * amount;
        DPV[acceptor.i][acceptor.j] += donor.v * amount;
    }
}

FaceState make_phase_state_from_cell(int i, int j,
                                     int Nx, int Ny, int fict_x, int fict_y,
                                     const Mask2D& solid_mask,
                                     const Field2& rho_tilde,
                                     const Field2& u_tilde,
                                     const Field2& v_tilde,
                                     const Field2& E) {
    FaceState state;
    if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return state;
    }
    state.rho = rho_tilde[i][j];
    state.u = u_tilde[i][j];
    state.v = v_tilde[i][j];
    state.E = E[i][j];
    state.blocked = false;
    state.physical = true;
    state.i = i;
    state.j = j;
    return state;
}

FaceState inflow_phase_state(const MaderConfig& config, double g) {
    FaceState state;
    const double rho = std::max(config.inflow_rho, EPS);
    const double I = std::max(config.inflow_p / (rho * (g - 1.0)), EPS);
    state.rho = rho;
    state.u = config.inflow_u;
    state.v = config.inflow_v;
    state.E = I + 0.5 * (state.u * state.u + state.v * state.v);
    state.blocked = false;
    state.physical = false;
    return state;
}

FaceState make_vertical_side(int i, int j,
                             bool is_left_side,
                             int Nx, int Ny, int fict_x, int fict_y,
                             const Mask2D& solid_mask,
                             const std::vector<unsigned char>& step_rows,
                             const Field2& rho_tilde,
                             const Field2& u_tilde,
                             const Field2& v_tilde,
                             const Field2& E,
                             const MaderConfig& config,
                             double g) {
    const int first_phys = fict_x;
    const int last_phys = Nx - fict_x - 1;
    const int cell_i = i;

    if (cell_i < first_phys) {
        if (step_rows[j]) {
            return FaceState{};
        }
        return inflow_phase_state(config, g);
    }
    if (cell_i > last_phys) {
        return make_phase_state_from_cell(last_phys, j, Nx, Ny, fict_x, fict_y,
                                          solid_mask, rho_tilde, u_tilde, v_tilde, E);
    }

    if (!is_fluid_cell(cell_i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return FaceState{};
    }

    return make_phase_state_from_cell(cell_i, j, Nx, Ny, fict_x, fict_y,
                                      solid_mask, rho_tilde, u_tilde, v_tilde, E);
}

FaceState make_horizontal_side(int i, int j,
                               int Nx, int Ny, int fict_x, int fict_y,
                               const Mask2D& solid_mask,
                               const Field2& rho_tilde,
                               const Field2& u_tilde,
                               const Field2& v_tilde,
                               const Field2& E) {
    const int first_phys = fict_y;
    const int last_phys = Ny - fict_y - 1;

    if (j < first_phys || j > last_phys) {
        return FaceState{};
    }
    if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
        return FaceState{};
    }
    return make_phase_state_from_cell(i, j, Nx, Ny, fict_x, fict_y,
                                      solid_mask, rho_tilde, u_tilde, v_tilde, E);
}

} // namespace

Mask2D build_step_mask(double x_min, double y_min, double dx, double dy,
                       int Nx, int Ny, int fict_x, int fict_y,
                       const MaderConfig& config) {
    Mask2D solid_mask(Nx, std::vector<unsigned char>(Ny, 0));
    if (!config.step_enabled) {
        return solid_mask;
    }

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        const double x = x_min + (i - fict_x + 0.5) * dx;
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            const double y = y_min + (j - fict_y + 0.5) * dy;
            if (x <= config.step_x_end + EPS && y <= config.step_y_end + EPS) {
                solid_mask[i][j] = 1;
            }
        }
    }
    return solid_mask;
}

void apply_step_mask(Field3& u,
                     const Mask2D& solid_mask,
                     const std::vector<double>& fill_cons) {
    for (size_t i = 0; i < u.size(); ++i) {
        for (size_t j = 0; j < u[i].size(); ++j) {
            if (solid_mask[i][j]) {
                u[i][j] = fill_cons;
            }
        }
    }
}

double calc_time_step_masked(const Field3& v_cons,
                             const Mask2D& solid_mask,
                             double dx, double dy, double cfl, double g,
                             int fict_x, int fict_y) {
    double max_sig_x = 0.0;
    double max_sig_y = 0.0;
    const int Nx = static_cast<int>(v_cons.size());
    const int Ny = static_cast<int>(v_cons[0].size());

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                continue;
            }
            const std::vector<double> prim = cons_to_noncons(v_cons[i][j], g);
            const double c = calc_sound_speed(prim, g);
            max_sig_x = std::max(max_sig_x, std::abs(prim[U]) + c);
            max_sig_y = std::max(max_sig_y, std::abs(prim[V]) + c);
        }
    }

    const double inv_dt = (max_sig_x / dx) + (max_sig_y / dy);
    if (inv_dt <= EPS) {
        return cfl * std::min(dx, dy);
    }
    return cfl / inv_dt;
}

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
                   const MaderConfig& config) {
    Field3 state = u_prev;
    const std::vector<double> inflow_cons = make_inflow_cons(config, g);
    const std::vector<unsigned char> step_rows = build_step_rows(solid_mask, Nx, Ny, fict_x, fict_y);

    apply_step_mask(state, solid_mask, inflow_cons);
    apply_mader_boundaries(state, Nx, Ny, fict_x, fict_y, solid_mask, config, g);

    Field2 rho = make_field2(Nx, Ny);
    Field2 u = make_field2(Nx, Ny);
    Field2 v = make_field2(Nx, Ny);
    Field2 p = make_field2(Nx, Ny);
    Field2 I = make_field2(Nx, Ny);

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                continue;
            }
            const std::vector<double> prim = cons_to_noncons(state[i][j], g);
            rho[i][j] = prim[RHO];
            u[i][j] = prim[U];
            v[i][j] = prim[V];
            p[i][j] = prim[P];
            I[i][j] = prim[P] / (prim[RHO] * (g - 1.0));
        }
    }

    Field2 qx = make_field2(Nx, Ny);
    Field2 qy = make_field2(Nx, Ny);

    for (int i = fict_x; i < Nx - fict_x - 1; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask) ||
                !is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }
            const double du = u[i][j] - u[i + 1][j];
            if (du > 0.0) {
                qx[i][j] = config.visc * 0.5 * (rho[i][j] + rho[i + 1][j]) * du;
            }
        }
    }

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y - 1; ++j) {
            if (!is_fluid_cell(i, j, Nx, Ny, fict_x, fict_y, solid_mask) ||
                !is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) {
                continue;
            }
            const double dv = v[i][j] - v[i][j + 1];
            if (dv > 0.0) {
                qy[i][j] = config.visc * 0.5 * (rho[i][j] + rho[i][j + 1]) * dv;
            }
        }
    }

    Field2 u_tilde = u;
    Field2 v_tilde = v;

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                continue;
            }

            const std::vector<double> left_cons = sample_cons(state, i, j, i - 1, j,
                                                              Nx, Ny, fict_x, fict_y,
                                                              solid_mask, step_rows, inflow_cons, g);
            const std::vector<double> right_cons = sample_cons(state, i, j, i + 1, j,
                                                               Nx, Ny, fict_x, fict_y,
                                                               solid_mask, step_rows, inflow_cons, g);
            const std::vector<double> down_cons = sample_cons(state, i, j, i, j - 1,
                                                              Nx, Ny, fict_x, fict_y,
                                                              solid_mask, step_rows, inflow_cons, g);
            const std::vector<double> up_cons = sample_cons(state, i, j, i, j + 1,
                                                            Nx, Ny, fict_x, fict_y,
                                                            solid_mask, step_rows, inflow_cons, g);

            const std::vector<double> left_prim = cons_to_noncons(left_cons, g);
            const std::vector<double> right_prim = cons_to_noncons(right_cons, g);
            const std::vector<double> down_prim = cons_to_noncons(down_cons, g);
            const std::vector<double> up_prim = cons_to_noncons(up_cons, g);

            const double q_left = (i > fict_x && is_fluid_cell(i - 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i - 1][j] : 0.0;
            const double q_right = (i < Nx - fict_x - 1 && is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i][j] : 0.0;
            const double q_down = (j > fict_y && is_fluid_cell(i, j - 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j - 1] : 0.0;
            const double q_up = (j < Ny - fict_y - 1 && is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j] : 0.0;

            u_tilde[i][j] = u[i][j] - (dt / std::max(rho[i][j], EPS) / dx) *
                                      ((right_prim[P] - left_prim[P]) + (q_right - q_left));
            v_tilde[i][j] = v[i][j] - (dt / std::max(rho[i][j], EPS) / dy) *
                                      ((up_prim[P] - down_prim[P]) + (q_up - q_down));
        }
    }

    Field2 rho_tilde = rho;
    Field2 I_tilde = I;
    Field2 E = make_field2(Nx, Ny);

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                continue;
            }

            const double U1 = face_u_sum_left(i, j, Nx, Ny, fict_x, fict_y, solid_mask, step_rows, u, u_tilde, config);
            const double U2 = face_u_sum_right(i, j, Nx, Ny, fict_x, fict_y, solid_mask, u, u_tilde);
            const double V1 = face_v_sum_down(i, j, Nx, Ny, fict_x, fict_y, solid_mask, v, v_tilde);
            const double V2 = face_v_sum_up(i, j, Nx, Ny, fict_x, fict_y, solid_mask, v, v_tilde);
            const double T3 = u[i][j] + u_tilde[i][j];
            const double T1 = v[i][j] + v_tilde[i][j];

            const double q_left = (i > fict_x && is_fluid_cell(i - 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i - 1][j] : 0.0;
            const double q_right = (i < Nx - fict_x - 1 && is_fluid_cell(i + 1, j, Nx, Ny, fict_x, fict_y, solid_mask)) ? qx[i][j] : 0.0;
            const double q_down = (j > fict_y && is_fluid_cell(i, j - 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j - 1] : 0.0;
            const double q_up = (j < Ny - fict_y - 1 && is_fluid_cell(i, j + 1, Nx, Ny, fict_x, fict_y, solid_mask)) ? qy[i][j] : 0.0;

            rho_tilde[i][j] = rho[i][j] - (rho[i][j] * dt * 0.5) *
                                        ((U2 - U1) / dx + (V2 - V1) / dy);
            rho_tilde[i][j] = std::max(rho_tilde[i][j], EPS);

            const double zip_term =
                (p[i][j] / dx) * (U2 - U1) +
                (q_right / dx) * (U2 - T3) +
                (q_left / dx) * (T3 - U1) +
                (p[i][j] / dy) * (V2 - V1) +
                (q_up / dy) * (V2 - T1) +
                (q_down / dy) * (T1 - V1);

            I_tilde[i][j] = I[i][j] - (dt / (4.0 * std::max(rho[i][j], EPS))) * zip_term;
            I_tilde[i][j] = std::max(I_tilde[i][j], EPS);
            E[i][j] = I_tilde[i][j] + 0.5 * (u_tilde[i][j] * u_tilde[i][j] + v_tilde[i][j] * v_tilde[i][j]);
        }
    }

    Field2 DM = make_field2(Nx, Ny);
    Field2 DE = make_field2(Nx, Ny);
    Field2 DPU = make_field2(Nx, Ny);
    Field2 DPV = make_field2(Nx, Ny);

    for (int iface = fict_x; iface <= Nx - fict_x; ++iface) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            const int il = iface - 1;
            const int ir = iface;
            FaceState left = make_vertical_side(il, j, true, Nx, Ny, fict_x, fict_y,
                                                solid_mask, step_rows, rho_tilde, u_tilde, v_tilde, E, config, g);
            FaceState right = make_vertical_side(ir, j, false, Nx, Ny, fict_x, fict_y,
                                                 solid_mask, step_rows, rho_tilde, u_tilde, v_tilde, E, config, g);
            if (left.blocked || right.blocked) {
                continue;
            }
            if (!left.physical && !right.physical) {
                continue;
            }

            double denom = 1.0 + (left.u - right.u) * dt / dx;
            if (std::abs(denom) < EPS) {
                denom = (denom >= 0.0) ? EPS : -EPS;
            }
            double alpha = 0.5 * (left.u + right.u) * dt / dx / denom;
            alpha = clamp_fraction(alpha);

            const FaceState& donor = (alpha >= 0.0) ? left : right;
            const FaceState& acceptor = (alpha >= 0.0) ? right : left;
            const double amount = donor.rho * std::abs(alpha);
            accumulate_transport(DM, DE, DPU, DPV, donor, acceptor, amount);
        }
    }

    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int jface = fict_y; jface <= Ny - fict_y; ++jface) {
            const int jd = jface - 1;
            const int ju = jface;
            FaceState down = make_horizontal_side(i, jd, Nx, Ny, fict_x, fict_y,
                                                  solid_mask, rho_tilde, u_tilde, v_tilde, E);
            FaceState up = make_horizontal_side(i, ju, Nx, Ny, fict_x, fict_y,
                                                solid_mask, rho_tilde, u_tilde, v_tilde, E);
            if (down.blocked || up.blocked) {
                continue;
            }

            double denom = 1.0 + (down.v - up.v) * dt / dy;
            if (std::abs(denom) < EPS) {
                denom = (denom >= 0.0) ? EPS : -EPS;
            }
            double beta = 0.5 * (down.v + up.v) * dt / dy / denom;
            beta = clamp_fraction(beta);

            const FaceState& donor = (beta >= 0.0) ? down : up;
            const FaceState& acceptor = (beta >= 0.0) ? up : down;
            const double amount = donor.rho * std::abs(beta);
            accumulate_transport(DM, DE, DPU, DPV, donor, acceptor, amount);
        }
    }

    u_next = state;
    for (int i = fict_x; i < Nx - fict_x; ++i) {
        for (int j = fict_y; j < Ny - fict_y; ++j) {
            if (solid_mask[i][j]) {
                u_next[i][j] = inflow_cons;
                continue;
            }

            const double rho_new = std::max(rho_tilde[i][j] + DM[i][j], EPS);
            const double mom_x = rho_tilde[i][j] * u_tilde[i][j] + DPU[i][j];
            const double mom_y = rho_tilde[i][j] * v_tilde[i][j] + DPV[i][j];
            const double u_new = mom_x / rho_new;
            const double v_new = mom_y / rho_new;
            const double E_new = (rho_tilde[i][j] * E[i][j] + DE[i][j]) / rho_new;
            const double I_new = std::max(E_new - 0.5 * (u_new * u_new + v_new * v_new), EPS);
            const double p_new = std::max((g - 1.0) * rho_new * I_new, EPS);

            u_next[i][j] = noncons_to_cons(std::vector<double>{rho_new, u_new, v_new, p_new}, g);
            enforce_physical_state(u_next[i][j], g);
        }
    }

    apply_step_mask(u_next, solid_mask, inflow_cons);
    apply_mader_boundaries(u_next, Nx, Ny, fict_x, fict_y, solid_mask, config, g);
}
