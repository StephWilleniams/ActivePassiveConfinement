#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <array>
#include <utility>
#include <yaml-cpp/yaml.h>

std::mt19937 gen;
std::uniform_real_distribution<> unif01(0.0, 1.0);
std::normal_distribution<> unif01n(0.0, 1.0);

enum class PotentialType { YUKAWA, WCA, ELASTIC, ENTRAINMENT };

constexpr double TWO_ONE_SIXTH = 1.122462048309373; // 2^(1/6) for WCA cut-off
constexpr double MIN_FORCE_DISTANCE = 1e-6;         // Distance floor to prevent division by zero
constexpr double ZERO_FORCE_THRESHOLD = 1e-12;      // Minimum force magnitude to process
const double SWIMMER_R_BASE = 2.0 / (4.0 + std::sqrt(2.0));
const double SWIMMER_BODY_LENGTH = 2.0;

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------
struct Config {
    // Simulation parameters
    int num_swimmers;
    int num_passives;
    int base_steps;
    double base_dt;
    double x_max;
    double y_max;
    double box_x;
    double inv_box_x;
    unsigned int seed;
    
    // Physics parameters
    double u0;
    double u_wake;
    double alpha;
    double alpha_rot;
    double alpha_rot_wall;
    double velocity;
    double rot_diffusion;
    double passive_diffusion;
    double passive_radius;
    double polarity;
    PotentialType potential;

    // Safety caps for numerical stability
    double force_cap;            // Ceiling on raw pairwise force (r->0 guard)
    double max_step_fraction;    // Steric displacement cap per step (fraction of active_diameter)
    double active_diameter;      // Reference lengthscale for max_step_fraction
    double max_angular_fraction; // Angular displacement cap per step (radians)

    // Boundary interaction parameters
    double boundary_threshold;
    double kick_rate;
    double kick_torque;

    // Precomputed drag equations
    double frot;
    double fpar;
    double fperp;
};

// -----------------------------------------------------------------------------
// Particle Data Structures
// -----------------------------------------------------------------------------
struct Swimmer {
    double x, y, theta, polarity;
    double Fx = 0.0, Fy = 0.0, torque = 0.0;
    double torque_wall = 0.0;
    
    double torque_kick = 0.0;
    double velocity, rot_diffusion;

    bool is_scattering = false;
    double scatter_torque_dir = 0.0;
    double scatter_elapsed = 0.0; 

    std::array<double, 3> radius;
    std::array<double, 3> segment_offsets; 

    Swimmer(double x_pos, double y_pos, double theta_pos, double pol_pos, double vel, double rot_diff) 
        : x(x_pos), y(y_pos), theta(theta_pos), polarity(pol_pos), velocity(vel), rot_diffusion(rot_diff) 
    {
        double r_base = SWIMMER_R_BASE;
        radius = { r_base, std::sqrt(2.0) * r_base, 2.0 * r_base };
        segment_offsets = {
            polarity * (1.0 - radius[0]), 
            polarity * (1.0 - 2.0 * radius[0]), 
            polarity * (1.0 - 2.0 * radius[0] - radius[1])
        };
    }
};

struct Passive {
    double x, y;
    double radius;
    double Fx = 0.0, Fy = 0.0;
    bool is_interacting = false;

    Passive(double x_pos, double y_pos, double r) : x(x_pos), y(y_pos), radius(r) {}
};

struct SegmentData {
    double x, y, dx_local, dy_local;
};

// -----------------------------------------------------------------------------
// Neighbor List (Verlet list with skin distance)
// -----------------------------------------------------------------------------
struct NeighborList {
    std::vector<std::pair<int, int>> ss_pairs; // swimmer-swimmer
    std::vector<std::pair<int, int>> sp_pairs; // swimmer-passive
    std::vector<std::pair<int, int>> pp_pairs; // passive-passive

    std::vector<double> ref_sx, ref_sy; // swimmer centers at last build
    std::vector<double> ref_px, ref_py; // passive centers at last build

    double skin = 0.0;
    double cutoff_ss_sq = 0.0;
    double cutoff_sp_sq = 0.0;
    double cutoff_pp_sq = 0.0;
    bool built = false;
};

// -----------------------------------------------------------------------------
// Utility Functions
// -----------------------------------------------------------------------------
double bounded_rand(double lower, double upper) {
    return lower + unif01(gen) * (upper - lower);
}

// Fast power of 6
inline double pow6(double x) {
    double x2 = x * x;
    return x2 * x2 * x2;
}

// Optimised Minimum Image Convention using precalculated inverse box size
inline double apply_pbc_dx(double dx, double box_x, double inv_box_x) {
    double half = 0.5 * box_x;
    if (dx >  half) return dx - box_x;
    if (dx < -half) return dx + box_x;
    return dx;
}

// Fast modulo-free Periodic Boundary Conditions and boundary clamping
inline void apply_pbc_and_clamp(double& x, double& y, const Config& cfg) {
    if (x > cfg.x_max) x -= cfg.box_x;
    else if (x < -cfg.x_max) x += cfg.box_x;

    if (y > cfg.y_max) y = cfg.y_max;
    else if (y < -cfg.y_max) y = -cfg.y_max;
}

// -----------------------------------------------------------------------------
// Physics Calculations
// -----------------------------------------------------------------------------

// Squared interaction cutoff per potential to avoid expensive sqrt() calls on far pairs
template <PotentialType PT>
inline double cutoff_sq(double sigma) {
    if constexpr (PT == PotentialType::WCA) { 
        double rc = TWO_ONE_SIXTH * sigma; 
        return rc * rc; 
    } else if constexpr (PT == PotentialType::ELASTIC) { 
        return sigma * sigma; 
    } else { 
        double rc = 3.0 * sigma; 
        return rc * rc; // YUKAWA, ENTRAINMENT
    }
}

// Templated force calculation to avoid branching in the inner loop
template <PotentialType PT>
inline double calculate_force_magnitude(double r, double sigma, const Config& cfg) {
    if (r < MIN_FORCE_DISTANCE) return cfg.force_cap; 
    
    double force_mag = 0.0;
    
    if constexpr (PT == PotentialType::YUKAWA) {
        double rho = r / sigma;
        if (rho < 3.0) { 
            double inv_rho = 1.0 / rho;
            force_mag = (cfg.u0 / sigma) * std::exp(-rho) * (inv_rho * inv_rho + 2.0 * inv_rho * inv_rho * inv_rho);
        }
    } else if constexpr (PT == PotentialType::WCA) {
        double r_cut = TWO_ONE_SIXTH * sigma;
        if (r < r_cut) {
            double sr6 = pow6(sigma / r);
            force_mag = 24.0 * cfg.u0 * sr6 * (2.0 * sr6 - 1.0) / r;
        }
    } else if constexpr (PT == PotentialType::ELASTIC) {
        if (r < sigma) {
            force_mag = cfg.u0 * (sigma - r);
        }
    } else if constexpr (PT == PotentialType::ENTRAINMENT) {
        double r_cut = TWO_ONE_SIXTH * sigma;
        if (r < r_cut) {
            double sr6 = pow6(sigma / r);
            force_mag += 24.0 * cfg.u0 * sr6 * (2.0 * sr6 - 1.0) / r;
        }
        if (r >= r_cut && r < 3.0 * sigma) {
            double rho = r / sigma;
            double inv_rho = 1.0 / rho;
            force_mag -= (cfg.u_wake / sigma) * std::exp(-rho) * (inv_rho * inv_rho + 2.0 * inv_rho * inv_rho * inv_rho);
        }
    }

    // Apply numerical guard capping the raw force
    if (force_mag > cfg.force_cap) {
        return cfg.force_cap;
    }

    return force_mag;
}

// (Re)builds the neighbour list if it's possible that any pair could have moved
// within interaction range since the last build. Safe for any potential type.
template <PotentialType PT>
void maybe_rebuild_neighbor_list(const std::vector<Swimmer>& swimmers, const std::vector<Passive>& passives, const Config& cfg, NeighborList& nl) {
    const size_t n_swimmers = swimmers.size();
    const size_t n_passives = passives.size();

    if (nl.built) {
        double max_disp = 0.0;
        for (size_t i = 0; i < n_swimmers; ++i) {
            double dx = apply_pbc_dx(swimmers[i].x - nl.ref_sx[i], cfg.box_x, cfg.inv_box_x);
            double dy = swimmers[i].y - nl.ref_sy[i];
            max_disp = std::max(max_disp, std::sqrt(dx * dx + dy * dy));
        }
        for (size_t i = 0; i < n_passives; ++i) {
            double dx = apply_pbc_dx(passives[i].x - nl.ref_px[i], cfg.box_x, cfg.inv_box_x);
            double dy = passives[i].y - nl.ref_py[i];
            max_disp = std::max(max_disp, std::sqrt(dx * dx + dy * dy));
        }
        if (2.0 * max_disp < nl.skin) return;
    }

    if (!nl.built) {
        double max_radius = (n_swimmers > 0) ? *std::max_element(swimmers[0].radius.begin(), swimmers[0].radius.end()) : 0.0;
        double max_seg_extent = 0.0;
        
        if (n_swimmers > 0) {
            for (double off : swimmers[0].segment_offsets) max_seg_extent = std::max(max_seg_extent, std::abs(off));
        }

        double range_ss = std::sqrt(cutoff_sq<PT>(2.0 * max_radius));
        double range_sp = std::sqrt(cutoff_sq<PT>(max_radius + cfg.passive_radius));
        double range_pp = std::sqrt(cutoff_sq<PT>(2.0 * cfg.passive_radius));

        nl.skin = 0.5 * cfg.active_diameter;
        double list_radius_ss = range_ss + 2.0 * max_seg_extent + nl.skin;
        double list_radius_sp = range_sp + max_seg_extent + nl.skin;
        double list_radius_pp = range_pp + nl.skin;

        nl.cutoff_ss_sq = list_radius_ss * list_radius_ss;
        nl.cutoff_sp_sq = list_radius_sp * list_radius_sp;
        nl.cutoff_pp_sq = list_radius_pp * list_radius_pp;

        nl.ref_sx.resize(n_swimmers); nl.ref_sy.resize(n_swimmers);
        nl.ref_px.resize(n_passives); nl.ref_py.resize(n_passives);
    }

    nl.ss_pairs.clear();
    for (size_t i = 0; i < n_swimmers; ++i) {
        for (size_t j = i + 1; j < n_swimmers; ++j) {
            double dx = apply_pbc_dx(swimmers[i].x - swimmers[j].x, cfg.box_x, cfg.inv_box_x);
            double dy = swimmers[i].y - swimmers[j].y;
            if (dx * dx + dy * dy < nl.cutoff_ss_sq) nl.ss_pairs.emplace_back((int)i, (int)j);
        }
    }

    nl.sp_pairs.clear();
    for (size_t i = 0; i < n_swimmers; ++i) {
        for (size_t j = 0; j < n_passives; ++j) {
            double dx = apply_pbc_dx(swimmers[i].x - passives[j].x, cfg.box_x, cfg.inv_box_x);
            double dy = swimmers[i].y - passives[j].y;
            if (dx * dx + dy * dy < nl.cutoff_sp_sq) nl.sp_pairs.emplace_back((int)i, (int)j);
        }
    }

    nl.pp_pairs.clear();
    for (size_t i = 0; i < n_passives; ++i) {
        for (size_t j = i + 1; j < n_passives; ++j) {
            double dx = apply_pbc_dx(passives[i].x - passives[j].x, cfg.box_x, cfg.inv_box_x);
            double dy = passives[i].y - passives[j].y;
            if (dx * dx + dy * dy < nl.cutoff_pp_sq) nl.pp_pairs.emplace_back((int)i, (int)j);
        }
    }

    for (size_t i = 0; i < n_swimmers; ++i) { nl.ref_sx[i] = swimmers[i].x; nl.ref_sy[i] = swimmers[i].y; }
    for (size_t i = 0; i < n_passives; ++i) { nl.ref_px[i] = passives[i].x; nl.ref_py[i] = passives[i].y; }
    nl.built = true;
}

template <PotentialType PT>
void calculate_potentials_impl(std::vector<Swimmer>& swimmers, std::vector<Passive>& passives, const Config& cfg, NeighborList& nl) {
    const size_t n_swimmers = swimmers.size();
    const size_t n_passives = passives.size();

    // Cache segment positions to avoid O(N^2) trigonometry evaluations
    static thread_local std::vector<std::array<SegmentData, 3>> segment_cache;
    segment_cache.resize(n_swimmers);

    // Initialise swimmer segments and reset forces
    for (size_t i = 0; i < n_swimmers; ++i) {
        swimmers[i].Fx = 0.0; swimmers[i].Fy = 0.0; swimmers[i].torque = 0.0;
        swimmers[i].torque_wall = 0.0;
        double ct = std::cos(swimmers[i].theta);
        double st = std::sin(swimmers[i].theta);

        for (int k = 0; k < 3; ++k) {
            segment_cache[i][k].dx_local = swimmers[i].segment_offsets[k] * ct;
            segment_cache[i][k].dy_local = swimmers[i].segment_offsets[k] * st;
            segment_cache[i][k].x = swimmers[i].x + segment_cache[i][k].dx_local;
            segment_cache[i][k].y = swimmers[i].y + segment_cache[i][k].dy_local;
        }
    }

    for (size_t i = 0; i < n_passives; ++i) {
        passives[i].Fx = 0.0; passives[i].Fy = 0.0;
        passives[i].is_interacting = false;
    }

    maybe_rebuild_neighbor_list<PT>(swimmers, passives, cfg, nl);

    // 1. Swimmer-Swimmer Interactions
    for (const auto& pr : nl.ss_pairs) {
        size_t i = (size_t)pr.first, j = (size_t)pr.second;
        #pragma GCC unroll 3
        for (int seg_i = 0; seg_i < 3; ++seg_i) {
            #pragma GCC unroll 3
            for (int seg_j = 0; seg_j < 3; ++seg_j) {

                double dx = apply_pbc_dx(segment_cache[i][seg_i].x - segment_cache[j][seg_j].x, cfg.box_x, cfg.inv_box_x);
                double dy = segment_cache[i][seg_i].y - segment_cache[j][seg_j].y;
                double r2 = dx * dx + dy * dy;
                double sigma = swimmers[i].radius[seg_i] + swimmers[j].radius[seg_j];

                if (r2 >= cutoff_sq<PT>(sigma)) continue;

                double r = std::sqrt(r2);
                double force_mag = calculate_force_magnitude<PT>(r, sigma, cfg);

                if (std::abs(force_mag) > ZERO_FORCE_THRESHOLD) {
                    double inv_r = 1.0 / r;
                    double fx = force_mag * dx * inv_r;
                    double fy = force_mag * dy * inv_r;

                    swimmers[i].Fx += fx;
                    swimmers[i].Fy += fy;
                    swimmers[i].torque += (segment_cache[i][seg_i].dx_local * fy - segment_cache[i][seg_i].dy_local * fx);

                    swimmers[j].Fx -= fx;
                    swimmers[j].Fy -= fy;
                    swimmers[j].torque -= (segment_cache[j][seg_j].dx_local * fy - segment_cache[j][seg_j].dy_local * fx);
                }
            }
        }
    }

    // 2. Swimmer-Passive Interactions 
    for (const auto& pr : nl.sp_pairs) {
        size_t i = (size_t)pr.first, j = (size_t)pr.second;
        #pragma GCC unroll 3
        for (int seg_i = 0; seg_i < 3; ++seg_i) {
            double dx = apply_pbc_dx(segment_cache[i][seg_i].x - passives[j].x, cfg.box_x, cfg.inv_box_x);
            double dy = segment_cache[i][seg_i].y - passives[j].y;
            double r2 = dx * dx + dy * dy;
            double sigma = swimmers[i].radius[seg_i] + passives[j].radius;

            if (r2 >= cutoff_sq<PT>(sigma)) continue;

            double r = std::sqrt(r2);
            double force_mag = calculate_force_magnitude<PT>(r, sigma, cfg);

            if (std::abs(force_mag) > ZERO_FORCE_THRESHOLD) {
                double inv_r = 1.0 / r;
                double fx = force_mag * dx * inv_r;
                double fy = force_mag * dy * inv_r;

                swimmers[i].Fx += fx;
                swimmers[i].Fy += fy;
                swimmers[i].torque += (segment_cache[i][seg_i].dx_local * fy - segment_cache[i][seg_i].dy_local * fx);

                passives[j].Fx -= fx;
                passives[j].Fy -= fy;
                passives[j].is_interacting = true;
            }
        }
    }

    // 3. Passive-Passive Interactions
    for (const auto& pr : nl.pp_pairs) {
        size_t i = (size_t)pr.first, j = (size_t)pr.second;
        
        double dx = apply_pbc_dx(passives[i].x - passives[j].x, cfg.box_x, cfg.inv_box_x);
        double dy = passives[i].y - passives[j].y;
        double r2 = dx * dx + dy * dy;
        double sigma = 2.0 * cfg.passive_radius; 

        if (r2 >= cutoff_sq<PT>(sigma)) continue;

        double r = std::sqrt(r2);
        double force_mag = calculate_force_magnitude<PT>(r, sigma, cfg);

        if (std::abs(force_mag) > ZERO_FORCE_THRESHOLD) {
            double inv_r = 1.0 / r;
            double fx = force_mag * dx * inv_r;
            double fy = force_mag * dy * inv_r;

            passives[i].Fx += fx;
            passives[i].Fy += fy;
            passives[j].Fx -= fx;
            passives[j].Fy -= fy;

            // passives[i].is_interacting = true;
            // passives[j].is_interacting = true;
        }
    }

    // 4. Y-Boundary Interactions
    for (size_t i = 0; i < n_swimmers; ++i) {
        #pragma GCC unroll 3
        for (int seg_i = 0; seg_i < 3; ++seg_i) {
            double sigma = swimmers[i].radius[seg_i];

            // Bottom Boundary
            double dist_bottom = segment_cache[i][seg_i].y + cfg.y_max;
            double f_bottom = (PT == PotentialType::ENTRAINMENT) 
                            ? calculate_force_magnitude<PotentialType::WCA>(dist_bottom, sigma, cfg)
                            : calculate_force_magnitude<PT>(dist_bottom, sigma, cfg);
            
            if (f_bottom > 0.0) {
                swimmers[i].Fy += f_bottom;
                swimmers[i].torque_wall += (segment_cache[i][seg_i].dx_local * f_bottom);
            }

            // Top Boundary
            double dist_top = cfg.y_max - segment_cache[i][seg_i].y;
            double f_top = (PT == PotentialType::ENTRAINMENT) 
                         ? calculate_force_magnitude<PotentialType::WCA>(dist_top, sigma, cfg)
                         : calculate_force_magnitude<PT>(dist_top, sigma, cfg);
            
            if (f_top > 0.0) {
                swimmers[i].Fy -= f_top;
                swimmers[i].torque_wall -= (segment_cache[i][seg_i].dx_local * f_top);
            }
        }
    }
    
    for (size_t i = 0; i < n_passives; ++i) {
        // Bottom Boundary
        double dist_bottom = passives[i].y + cfg.y_max;
        double f_bottom = (PT == PotentialType::ENTRAINMENT)
                        ? calculate_force_magnitude<PotentialType::WCA>(dist_bottom, passives[i].radius, cfg)
                        : calculate_force_magnitude<PT>(dist_bottom, passives[i].radius, cfg);
        
        if (f_bottom > 0.0) {
            passives[i].Fy += f_bottom;
        }
        
        // Top Boundary
        double dist_top = cfg.y_max - passives[i].y;
        double f_top = (PT == PotentialType::ENTRAINMENT)
                     ? calculate_force_magnitude<PotentialType::WCA>(dist_top, passives[i].radius, cfg)
                     : calculate_force_magnitude<PT>(dist_top, passives[i].radius, cfg);
        
        if (f_top > 0.0) {
            passives[i].Fy -= f_top;
        }
    }
}

// Wrapper to dispatch the template based on configuration
void calculate_potentials(std::vector<Swimmer>& swimmers, std::vector<Passive>& passives, const Config& cfg, NeighborList& nl) {
    switch (cfg.potential) {
        case PotentialType::YUKAWA:      calculate_potentials_impl<PotentialType::YUKAWA>(swimmers, passives, cfg, nl); break;
        case PotentialType::WCA:         calculate_potentials_impl<PotentialType::WCA>(swimmers, passives, cfg, nl); break;
        case PotentialType::ELASTIC:     calculate_potentials_impl<PotentialType::ELASTIC>(swimmers, passives, cfg, nl); break;
        case PotentialType::ENTRAINMENT: calculate_potentials_impl<PotentialType::ENTRAINMENT>(swimmers, passives, cfg, nl); break;
    }
}

// -----------------------------------------------------------------------------
// Movement & Kinematics
// -----------------------------------------------------------------------------

void apply_boundary_kicks(std::vector<Swimmer>& swimmers, double dt, const Config& cfg) {
    double prob = 1.0 - std::exp(-cfg.kick_rate * dt);

    for (auto& s : swimmers) {
        bool near_top = (cfg.y_max - s.y) <= cfg.boundary_threshold;
        bool near_bottom = (s.y - (-cfg.y_max)) <= cfg.boundary_threshold;

        if (!s.is_scattering) {
            // Trigger a kick if in the boundary zone and the Poisson event occurs
            if ((near_top || near_bottom) && unif01(gen) < prob) {
                s.is_scattering = true;
                s.scatter_elapsed = 0.0;
                double ct = std::cos(s.theta);

                if (near_top) {
                    // Left-moving (ct < 0) rotates CCW (+1), Right-moving (ct > 0) rotates CW (-1)
                    s.scatter_torque_dir = (ct < 0.0) ? 1.0 : -1.0;
                } else {
                    // Left-moving (ct < 0) rotates CW (-1), Right-moving (ct > 0) rotates CCW (+1)
                    s.scatter_torque_dir = (ct < 0.0) ? -1.0 : 1.0;
                }
            }
        } else {
            bool far_from_wall = (!near_top && !near_bottom);

            s.scatter_elapsed += dt;
            bool timed_out = s.scatter_elapsed >= (1.0 / s.rot_diffusion);
            double ct_now = std::cos(s.theta);
            double dir_now = 0.0;
            if (near_top)         dir_now = (ct_now < 0.0) ?  1.0 : -1.0;
            else if (near_bottom) dir_now = (ct_now < 0.0) ? -1.0 :  1.0;
            bool sign_flipped = (dir_now != 0.0) && (dir_now != s.scatter_torque_dir);

            if (far_from_wall || timed_out || sign_flipped) {
                s.is_scattering = false;
                s.scatter_torque_dir = 0.0;
                s.scatter_elapsed = 0.0;
            }
        }

        // Apply constant torque to the dedicated torque_kick term
        s.torque_kick = s.is_scattering ? cfg.kick_torque * s.scatter_torque_dir : 0.0;
    }
}

void move_particles(std::vector<Swimmer>& swimmers, std::vector<Passive>& passives, double dt, const Config& cfg) {
    double sq_2_rot_dt = std::sqrt(2.0 * cfg.rot_diffusion * dt);
    double sq_2_pas_dt = std::sqrt(2.0 * cfg.passive_diffusion * dt);
    double inv_fpar = 1.0 / cfg.fpar;
    double inv_fperp = 1.0 / cfg.fperp;

    // Determine absolute displacement bounds to prevent excessive unphysical leaps
    double max_speed = (cfg.max_step_fraction * cfg.active_diameter) / dt;
    double max_angular_speed = cfg.max_angular_fraction / dt;

    for (auto& s : swimmers) {
        double ct = std::cos(s.theta);
        double st = std::sin(s.theta);

        double f_par_local =  s.Fx * ct + s.Fy * st;
        double f_perp_local = -s.Fx * st + s.Fy * ct;

        double v_par_local = f_par_local * inv_fpar;
        double v_perp_local = f_perp_local * inv_fperp;

        // Ensure steric translation complies with the single-step speed cap
        double speed = std::sqrt(v_par_local * v_par_local + v_perp_local * v_perp_local);
        if (speed > max_speed) {
            double scale = max_speed / speed;
            v_par_local *= scale;
            v_perp_local *= scale;
        }

        double vx_steric = cfg.alpha * (v_par_local * ct - v_perp_local * st);
        double vy_steric = cfg.alpha * (v_par_local * st + v_perp_local * ct);

        s.x += (s.velocity * ct + vx_steric) * dt;
        s.y += (s.velocity * st + vy_steric) * dt;
        
        // Ensure steric angular change complies with the rotational displacement cap
        double omega_steric = (cfg.alpha_rot * s.torque
                               + cfg.alpha_rot_wall * s.torque_wall) / cfg.frot;
        if (std::abs(omega_steric) > max_angular_speed) {
            omega_steric = std::copysign(max_angular_speed, omega_steric);
        }
        
        // Combine capped steric mechanics with explicit scattering and random diffusion
        double omega = omega_steric + cfg.alpha_rot * s.torque_kick;
        s.theta += (omega * dt) + sq_2_rot_dt * unif01n(gen);

        // Angle normalisation 
        if (s.theta > 2.0 * M_PI) s.theta -= 2.0 * M_PI;
        else if (s.theta < 0.0) s.theta += 2.0 * M_PI;

        apply_pbc_and_clamp(s.x, s.y, cfg);
    }

    for (auto& p : passives) {
        double vx = cfg.alpha * p.Fx;
        double vy = cfg.alpha * p.Fy;

        // Ensure passive particle translation complies with speed caps
        double speed = std::sqrt(vx * vx + vy * vy);
        if (speed > max_speed) {
            double scale = max_speed / speed;
            vx *= scale;
            vy *= scale;
        }

        p.x += vx * dt + (sq_2_pas_dt * unif01n(gen));
        p.y += vy * dt + (sq_2_pas_dt * unif01n(gen));

        apply_pbc_and_clamp(p.x, p.y, cfg);
    }
}

// -----------------------------------------------------------------------------
// Output Generation
// -----------------------------------------------------------------------------
// Formats each row with snprintf into a growable buffer instead of iostream's
// locale-aware operator<<, which is dominated by per-call facet lookups and
// printf-style double formatting under the hood anyway. Produces byte-identical
// text (same "%.6f" rounding as std::fixed/setprecision(6)) at lower overhead.
void output_positions(const std::vector<Swimmer>& swimmers, const std::vector<Passive>& passives, int step, double dt, std::vector<char>& out_buf) {
    double time = step * dt;
    char line[160];

    for (const auto& s : swimmers) {
        int len = std::snprintf(line, sizeof(line), "%.6f, S, %.6f, %.6f, %.6f, 0\n", time, s.x, s.y, s.theta);
        out_buf.insert(out_buf.end(), line, line + len);
    }
    for (const auto& p : passives) {
        int len = std::snprintf(line, sizeof(line), "%.6f, P, %.6f, %.6f, 0.0, %d\n", time, p.x, p.y, p.is_interacting ? 1 : 0);
        out_buf.insert(out_buf.end(), line, line + len);
    }
}

// -----------------------------------------------------------------------------
// Main Execution
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // CLI: <multiplier> <record_stride> <config_path> <output_path>
    int multiplier = (argc < 2) ? 1 : std::stoi(argv[1]);
    int record_stride = (argc < 3) ? 1 : std::stoi(argv[2]);
    std::string config_path = (argc < 4) ? "c_code/main_code/config.yml" : argv[3];
    std::string output_path = (argc < 5) ? ("o_outputs/output_" + std::to_string(multiplier) + ".txt") : argv[4];

    Config cfg;
    try {
        YAML::Node config = YAML::LoadFile(config_path);
        
        cfg.num_swimmers       = config["simulation"]["num_swimmers"].as<int>();
        cfg.num_passives       = config["simulation"]["num_passives"].as<int>();
        cfg.base_steps         = config["simulation"]["base_steps"].as<int>();
        cfg.base_dt            = config["simulation"]["base_dt"].as<double>();
        cfg.x_max              = config["simulation"]["x_max"].as<double>();
        cfg.y_max              = config["simulation"]["y_max"].as<double>();
        cfg.seed               = config["simulation"]["seed"].as<unsigned int>();
        
        cfg.box_x              = 2.0 * cfg.x_max;
        cfg.inv_box_x          = 1.0 / cfg.box_x;

        cfg.u0                   = config["physics"]["u0"].as<double>();
        cfg.u_wake               = config["physics"]["u_wake"].as<double>();
        cfg.alpha                = config["physics"]["alpha"].as<double>();
        cfg.alpha_rot            = config["physics"]["alpha_rot"].as<double>();
        cfg.alpha_rot_wall       = config["physics"]["alpha_rot_wall"]
                                 ? config["physics"]["alpha_rot_wall"].as<double>()
                                 : cfg.alpha_rot;
        cfg.velocity             = config["physics"]["velocity"].as<double>();
        cfg.rot_diffusion        = config["physics"]["rot_diffusion"].as<double>();
        cfg.force_cap            = config["physics"]["force_cap"].as<double>();
        cfg.max_step_fraction    = config["physics"]["max_step_fraction"].as<double>();
        cfg.max_angular_fraction = config["physics"]["max_angular_fraction"].as<double>();
        cfg.passive_diffusion    = config["physics"]["passive_diffusion"].as<double>();
        cfg.passive_radius       = config["physics"]["passive_radius"].as<double>();
        cfg.polarity             = config["physics"]["polarity"].as<double>();

        // Reference lengthscale, derived from the same SWIMMER_R_BASE as the swimmer constructor
        cfg.active_diameter    = 2.0 * SWIMMER_R_BASE;

        cfg.boundary_threshold = config["physics"]["boundary_threshold"].as<double>();
        cfg.kick_rate          = config["physics"]["kick_rate"].as<double>();
        cfg.kick_torque        = config["physics"]["kick_torque"].as<double>();

        std::string pot_str = config["physics"]["potential"].as<std::string>();
        if (pot_str == "WCA") cfg.potential = PotentialType::WCA;
        else if (pot_str == "ELASTIC") cfg.potential = PotentialType::ELASTIC;
        else if (pot_str == "ENTRAINMENT") cfg.potential = PotentialType::ENTRAINMENT;
        else cfg.potential = PotentialType::YUKAWA;
        double aspRatio = 1.0 / SWIMMER_R_BASE;
        cfg.fpar = SWIMMER_BODY_LENGTH * 2.0 * M_PI / (std::log(aspRatio) - 0.207 + (0.980 / aspRatio) - (0.133 / (aspRatio * aspRatio)));
        cfg.fperp = SWIMMER_BODY_LENGTH * 4.0 * M_PI / (std::log(aspRatio) + 0.839 + (0.185 / aspRatio) + (0.233 / (aspRatio * aspRatio)));
        double dr_rot = std::log(aspRatio) - 0.662 + (0.917 / aspRatio) - (0.050 / (aspRatio * aspRatio));
        cfg.frot = M_PI * SWIMMER_BODY_LENGTH * SWIMMER_BODY_LENGTH * SWIMMER_BODY_LENGTH / (3.0 * dr_rot);



    } catch (const YAML::Exception& e) {
        std::cerr << "Error reading " << config_path << ": " << e.what() << "\n";
        return 1;
    }

    gen.seed(cfg.seed);

    int N_steps = cfg.base_steps * multiplier;
    double dt = cfg.base_dt / multiplier;

    FILE* myfile = std::fopen(output_path.c_str(), "w");

    if (!myfile) {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }
    std::vector<char> out_buf;
    out_buf.reserve(1 << 20);
    constexpr size_t OUT_BUF_FLUSH_THRESHOLD = 1 << 20;
    auto flush_out_buf = [&]() {
        if (!out_buf.empty()) {
            std::fwrite(out_buf.data(), 1, out_buf.size(), myfile);
            out_buf.clear();
        }
    };

    std::vector<Swimmer> swimmers;
    swimmers.reserve(cfg.num_swimmers);
    for (int i = 0; i < cfg.num_swimmers; ++i) {
        swimmers.emplace_back(
            bounded_rand(-cfg.x_max, cfg.x_max), 
            bounded_rand(-cfg.y_max, cfg.y_max), 
            bounded_rand(0.0, 2.0 * M_PI),
            cfg.polarity,
            cfg.velocity,
            cfg.rot_diffusion
        );
    }
    
    std::vector<Passive> passives;
    passives.reserve(cfg.num_passives);
    for (int i = 0; i < cfg.num_passives; ++i) {
        passives.emplace_back(
            bounded_rand(-cfg.x_max, cfg.x_max), 
            bounded_rand(-cfg.y_max + cfg.passive_radius, cfg.y_max - cfg.passive_radius), 
            cfg.passive_radius
        );
    }

    output_positions(swimmers, passives, 0, dt, out_buf);

    NeighborList nl;

    // Main integration loop
    for (int step = 0; step < N_steps; ++step) {
        calculate_potentials(swimmers, passives, cfg, nl);
        apply_boundary_kicks(swimmers, dt, cfg);
        move_particles(swimmers, passives, dt, cfg);

        if ((step + 1) % multiplier == 0) {
            int frame = (step + 1) / multiplier;
            if (frame % record_stride == 0) {
                output_positions(swimmers, passives, step + 1, dt, out_buf);
                if (out_buf.size() >= OUT_BUF_FLUSH_THRESHOLD) flush_out_buf();
            }
        }
    }

    flush_out_buf();
    std::fclose(myfile);

    return 0;
}
