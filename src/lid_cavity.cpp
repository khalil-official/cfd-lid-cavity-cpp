#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
constexpr double PI = 3.141592653589793238462643383279502884;

struct Field {
    int rows{};
    int cols{};
    std::vector<double> values;

    Field() = default;
    Field(int r, int c, double value = 0.0)
        : rows(r), cols(c), values(static_cast<std::size_t>(r) * static_cast<std::size_t>(c), value) {}

    double& operator()(int i, int j) {
        return values[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols) + static_cast<std::size_t>(j)];
    }
    double operator()(int i, int j) const {
        return values[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols) + static_cast<std::size_t>(j)];
    }
};

struct Config {
    double lid_velocity = 1.0;
    double length = 1.0;

    int max_iterations = 30000;
    int minimum_iterations = 200;
    int consecutive_passes = 20;
    int maximum_pressure_failures = 3;

    double velocity_tolerance = 1e-8;
    double divergence_linf_tolerance = 1e-9;
    double divergence_l2_tolerance = 2e-10;
    double diverged_limit = 1e6;

    double cfl = 0.60;
    double dt_max = 0.01;
    double dt_min = 1e-8;
    double momentum_relaxation = 0.90;
    double pressure_relaxation = 1.0;

    int poisson_max_iterations = 5000;
    int poisson_check_every = 20;
    double poisson_absolute_tolerance = 1e-10;
    double poisson_relative_tolerance = 1e-9;
    double sor_omega = 0.0;

    int stagnation_window = 1500;
    double stagnation_minimum_reduction = 0.005;

    bool save_fields = true;
    bool strict_exit = false;
    std::string data_directory = "results/data";
};

struct PoissonInfo {
    int iterations = 0;
    bool converged = false;
    double absolute_residual = std::numeric_limits<double>::infinity();
    double relative_residual = std::numeric_limits<double>::infinity();
    double omega = 1.0;
};

struct Residuals {
    double velocity_update_linf = 0.0;
    double divergence_linf = 0.0;
    double divergence_l2 = 0.0;
    double global_mass_imbalance = 0.0;
};

struct GhiaMetrics {
    bool available = false;
    bool passed = false;
    double u_l2 = std::numeric_limits<double>::quiet_NaN();
    double v_l2 = std::numeric_limits<double>::quiet_NaN();
    double u_linf = std::numeric_limits<double>::quiet_NaN();
    double v_linf = std::numeric_limits<double>::quiet_NaN();
    double u_limit = std::numeric_limits<double>::quiet_NaN();
    double v_limit = std::numeric_limits<double>::quiet_NaN();
};

struct Result {
    int case_id = 0;
    int cells = 0;
    int reynolds = 0;
    std::string scheme;
    std::string pressure_solver;
    std::string status = "max_iterations";
    std::string quality = "needs_improvement";

    int iterations = 0;
    int local_max_iterations = 0;
    int failed_pressure_solves = 0;
    int consecutive_pass_count = 0;
    double runtime_seconds = 0.0;

    Field u_face;
    Field v_face;
    Field pressure;
    Field u_center;
    Field v_center;
    Field speed;
    Field vorticity;
    std::vector<double> x;
    std::vector<double> y;

    std::vector<double> velocity_history;
    std::vector<double> divergence_linf_history;
    std::vector<double> divergence_l2_history;
    std::vector<double> mass_history;
    std::vector<double> dt_history;
    std::vector<double> poisson_relative_history;
    std::vector<int> poisson_iteration_history;
    std::vector<bool> poisson_converged_history;

    Residuals final_residuals;
    double average_poisson_iterations = 0.0;
    double average_poisson_relative_residual = 0.0;
    double pressure_saturation_ratio = 0.0;
    GhiaMetrics ghia;
};

struct InitialState {
    bool available = false;
    Field u_face;
    Field v_face;
    Field pressure;
};

struct GhiaData {
    std::vector<double> y_u;
    std::vector<double> u;
    std::vector<double> x_v;
    std::vector<double> v;
};

static std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

static std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

static double maximum_absolute(const Field& field) {
    double result = 0.0;
    for (const double value : field.values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

static bool all_finite(const Field& field) {
    return std::all_of(field.values.begin(), field.values.end(), [](double value) {
        return std::isfinite(value);
    });
}

static void remove_mean(Field& field) {
    if (field.values.empty()) {
        return;
    }
    const double mean = std::accumulate(field.values.begin(), field.values.end(), 0.0)
                      / static_cast<double>(field.values.size());
    for (double& value : field.values) {
        value -= mean;
    }
}

static double u_with_wall_ghost(const Field& u, int i, int j, int cells, double lid_velocity) {
    if (i < 0) {
        return -u(0, j);
    }
    if (i >= cells) {
        return 2.0 * lid_velocity - u(cells - 1, j);
    }
    return u(i, j);
}

static double v_with_wall_ghost(const Field& v, int i, int j, int cells) {
    if (j < 0) {
        return -v(i, 0);
    }
    if (j >= cells) {
        return -v(i, cells - 1);
    }
    return v(i, j);
}

static Field divergence(const Field& u, const Field& v, double dx, double dy, int cells) {
    Field result(cells, cells, 0.0);
    for (int i = 0; i < cells; ++i) {
        for (int j = 0; j < cells; ++j) {
            result(i, j) = (u(i, j + 1) - u(i, j)) / dx
                         + (v(i + 1, j) - v(i, j)) / dy;
        }
    }
    return result;
}

static double poisson_residual(const Field& pressure, const Field& rhs, double h) {
    const int cells = pressure.rows;
    double maximum = 0.0;
    for (int i = 0; i < cells; ++i) {
        for (int j = 0; j < cells; ++j) {
            double laplacian = 0.0;
            if (i > 0) {
                laplacian += pressure(i - 1, j) - pressure(i, j);
            }
            if (i + 1 < cells) {
                laplacian += pressure(i + 1, j) - pressure(i, j);
            }
            if (j > 0) {
                laplacian += pressure(i, j - 1) - pressure(i, j);
            }
            if (j + 1 < cells) {
                laplacian += pressure(i, j + 1) - pressure(i, j);
            }
            laplacian /= h * h;
            maximum = std::max(maximum, std::abs(laplacian - rhs(i, j)));
        }
    }
    return maximum;
}

static std::pair<Field, PoissonInfo> solve_pressure_poisson(
    Field rhs,
    double h,
    const std::string& solver_name,
    const Config& config
) {
    const int cells = rhs.rows;
    const std::string method = upper(solver_name);
    const bool use_sor = method == "RBSOR";
    if (!use_sor && method != "RBGS") {
        throw std::runtime_error("Pressure solver must be RBGS or RBSOR");
    }

    const double rhs_mean = std::accumulate(rhs.values.begin(), rhs.values.end(), 0.0)
                          / static_cast<double>(rhs.values.size());
    double rhs_norm = 0.0;
    for (double& value : rhs.values) {
        value -= rhs_mean;
        rhs_norm = std::max(rhs_norm, std::abs(value));
    }
    rhs_norm = std::max(rhs_norm, 1e-30);

    Field pressure(cells, cells, 0.0);
    double omega = use_sor ? config.sor_omega : 1.0;
    if (use_sor && omega <= 0.0) {
        omega = 2.0 / (1.0 + std::sin(PI / static_cast<double>(cells)));
        omega = std::clamp(omega, 1.0, 1.95);
    }

    PoissonInfo info;
    info.omega = omega;

    for (int iteration = 1; iteration <= config.poisson_max_iterations; ++iteration) {
        for (int color = 0; color < 2; ++color) {
            for (int i = 0; i < cells; ++i) {
                // Visit only this color, preserving the previous update order.
                for (int j = (color + i) & 1; j < cells; j += 2) {
                    double neighbor_sum = 0.0;
                    int neighbor_count = 0;
                    if (i > 0) {
                        neighbor_sum += pressure(i - 1, j);
                        ++neighbor_count;
                    }
                    if (i + 1 < cells) {
                        neighbor_sum += pressure(i + 1, j);
                        ++neighbor_count;
                    }
                    if (j > 0) {
                        neighbor_sum += pressure(i, j - 1);
                        ++neighbor_count;
                    }
                    if (j + 1 < cells) {
                        neighbor_sum += pressure(i, j + 1);
                        ++neighbor_count;
                    }

                    const double candidate = (neighbor_sum - rhs(i, j) * h * h)
                                           / static_cast<double>(neighbor_count);
                    pressure(i, j) = use_sor
                        ? (1.0 - omega) * pressure(i, j) + omega * candidate
                        : candidate;
                }
            }
        }

        if (iteration % 50 == 0) {
            remove_mean(pressure);
        }

        if (iteration == 1
            || iteration % config.poisson_check_every == 0
            || iteration == config.poisson_max_iterations) {
            const double absolute = poisson_residual(pressure, rhs, h);
            const double relative = absolute / rhs_norm;
            info.iterations = iteration;
            info.absolute_residual = absolute;
            info.relative_residual = relative;
            if (absolute <= config.poisson_absolute_tolerance
                || relative <= config.poisson_relative_tolerance) {
                info.converged = true;
                remove_mean(pressure);
                return {pressure, info};
            }
        }
    }

    remove_mean(pressure);
    return {pressure, info};
}

static double compute_time_step(
    const Field& u,
    const Field& v,
    double h,
    double viscosity,
    const Config& config
) {
    const double maximum_velocity = std::max({
        maximum_absolute(u),
        maximum_absolute(v),
        config.lid_velocity,
        1e-12
    });
    const double convection_limit = config.cfl * h / maximum_velocity;
    const double diffusion_limit = 0.24 * h * h / std::max(viscosity, 1e-30);
    return std::clamp(
        std::min({convection_limit, diffusion_limit, config.dt_max}),
        config.dt_min,
        config.dt_max
    );
}

static std::tuple<Field, Field, double> predict_velocity(
    const Field& u,
    const Field& v,
    const Field& pressure,
    int reynolds,
    const std::string& scheme_name,
    const Config& config
) {
    const int cells = pressure.rows;
    const double h = config.length / static_cast<double>(cells);
    const double viscosity = config.lid_velocity * config.length / static_cast<double>(reynolds);
    const double dt = compute_time_step(u, v, h, viscosity, config);
    const std::string scheme = lower(scheme_name);

    Field u_star = u;
    Field v_star = v;

    for (int i = 0; i < cells; ++i) {
        for (int j = 1; j < cells; ++j) {
            const double u_center = u(i, j);
            const double v_at_u = 0.25 * (
                v(i, j - 1) + v(i + 1, j - 1)
                + v(i, j) + v(i + 1, j)
            );

            const double u_west = u(i, j - 1);
            const double u_east = u(i, j + 1);
            const double u_south = u_with_wall_ghost(u, i - 1, j, cells, config.lid_velocity);
            const double u_north = u_with_wall_ghost(u, i + 1, j, cells, config.lid_velocity);

            double du_dx = 0.0;
            double du_dy = 0.0;
            if (scheme == "central") {
                du_dx = (u_east - u_west) / (2.0 * h);
                du_dy = (u_north - u_south) / (2.0 * h);
            } else if (scheme == "upwind") {
                du_dx = u_center >= 0.0
                    ? (u_center - u_west) / h
                    : (u_east - u_center) / h;
                du_dy = v_at_u >= 0.0
                    ? (u_center - u_south) / h
                    : (u_north - u_center) / h;
            } else {
                throw std::runtime_error("Convection scheme must be upwind or central");
            }

            const double laplacian = (
                u_east - 2.0 * u_center + u_west
                + u_north - 2.0 * u_center + u_south
            ) / (h * h);
            const double pressure_gradient = (pressure(i, j) - pressure(i, j - 1)) / h;

            u_star(i, j) = u_center + config.momentum_relaxation * dt * (
                -u_center * du_dx - v_at_u * du_dy
                - pressure_gradient + viscosity * laplacian
            );
        }
    }

    for (int i = 1; i < cells; ++i) {
        for (int j = 0; j < cells; ++j) {
            const double v_center = v(i, j);
            const double u_at_v = 0.25 * (
                u(i - 1, j) + u(i - 1, j + 1)
                + u(i, j) + u(i, j + 1)
            );

            const double v_south = v(i - 1, j);
            const double v_north = v(i + 1, j);
            const double v_west = v_with_wall_ghost(v, i, j - 1, cells);
            const double v_east = v_with_wall_ghost(v, i, j + 1, cells);

            double dv_dx = 0.0;
            double dv_dy = 0.0;
            if (scheme == "central") {
                dv_dx = (v_east - v_west) / (2.0 * h);
                dv_dy = (v_north - v_south) / (2.0 * h);
            } else if (scheme == "upwind") {
                dv_dx = u_at_v >= 0.0
                    ? (v_center - v_west) / h
                    : (v_east - v_center) / h;
                dv_dy = v_center >= 0.0
                    ? (v_center - v_south) / h
                    : (v_north - v_center) / h;
            }

            const double laplacian = (
                v_east - 2.0 * v_center + v_west
                + v_north - 2.0 * v_center + v_south
            ) / (h * h);
            const double pressure_gradient = (pressure(i, j) - pressure(i - 1, j)) / h;

            v_star(i, j) = v_center + config.momentum_relaxation * dt * (
                -u_at_v * dv_dx - v_center * dv_dy
                - pressure_gradient + viscosity * laplacian
            );
        }
    }

    return {u_star, v_star, dt};
}

static Residuals calculate_residuals(
    const Field& u,
    const Field& v,
    const Field& old_u,
    const Field& old_v,
    double h,
    const Config& config
) {
    Residuals residuals;
    for (std::size_t k = 0; k < u.values.size(); ++k) {
        residuals.velocity_update_linf = std::max(
            residuals.velocity_update_linf,
            std::abs(u.values[k] - old_u.values[k]) / config.lid_velocity
        );
    }
    for (std::size_t k = 0; k < v.values.size(); ++k) {
        residuals.velocity_update_linf = std::max(
            residuals.velocity_update_linf,
            std::abs(v.values[k] - old_v.values[k]) / config.lid_velocity
        );
    }

    const Field div = divergence(u, v, h, h, u.rows);
    double square_sum = 0.0;
    for (const double value : div.values) {
        const double dimensionless = value / (config.lid_velocity / config.length);
        residuals.divergence_linf = std::max(residuals.divergence_linf, std::abs(dimensionless));
        square_sum += dimensionless * dimensionless;
    }
    residuals.divergence_l2 = std::sqrt(square_sum / static_cast<double>(div.values.size()));

    double boundary_flux = 0.0;
    const int cells = u.rows;
    for (int i = 0; i < cells; ++i) {
        boundary_flux += (u(i, cells) - u(i, 0)) * h;
    }
    for (int j = 0; j < cells; ++j) {
        boundary_flux += (v(cells, j) - v(0, j)) * h;
    }
    residuals.global_mass_imbalance = std::abs(boundary_flux)
                                    / (config.lid_velocity * config.length);
    return residuals;
}

static void build_cell_center_fields(Result& result, const Config& config) {
    const int cells = result.cells;
    const double h = config.length / static_cast<double>(cells);
    result.u_center = Field(cells, cells, 0.0);
    result.v_center = Field(cells, cells, 0.0);
    result.speed = Field(cells, cells, 0.0);
    result.vorticity = Field(cells, cells, 0.0);
    result.x.resize(static_cast<std::size_t>(cells));
    result.y.resize(static_cast<std::size_t>(cells));

    for (int k = 0; k < cells; ++k) {
        result.x[static_cast<std::size_t>(k)] = (static_cast<double>(k) + 0.5) * h;
        result.y[static_cast<std::size_t>(k)] = (static_cast<double>(k) + 0.5) * h;
    }

    for (int i = 0; i < cells; ++i) {
        for (int j = 0; j < cells; ++j) {
            result.u_center(i, j) = 0.5 * (result.u_face(i, j) + result.u_face(i, j + 1));
            result.v_center(i, j) = 0.5 * (result.v_face(i, j) + result.v_face(i + 1, j));
            result.speed(i, j) = std::hypot(result.u_center(i, j), result.v_center(i, j));
        }
    }

    for (int i = 1; i + 1 < cells; ++i) {
        for (int j = 1; j + 1 < cells; ++j) {
            const double dv_dx = (result.v_center(i, j + 1) - result.v_center(i, j - 1)) / (2.0 * h);
            const double du_dy = (result.u_center(i + 1, j) - result.u_center(i - 1, j)) / (2.0 * h);
            result.vorticity(i, j) = dv_dx - du_dy;
        }
    }
}

static double interpolate(const std::vector<double>& coordinates, const std::vector<double>& values, double query) {
    if (coordinates.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (query <= coordinates.front()) {
        return values.front();
    }
    if (query >= coordinates.back()) {
        return values.back();
    }
    const auto upper_it = std::upper_bound(coordinates.begin(), coordinates.end(), query);
    const std::size_t upper_index = static_cast<std::size_t>(std::distance(coordinates.begin(), upper_it));
    const std::size_t lower_index = upper_index - 1;
    const double fraction = (query - coordinates[lower_index])
                          / (coordinates[upper_index] - coordinates[lower_index]);
    return values[lower_index] + fraction * (values[upper_index] - values[lower_index]);
}

static GhiaData ghia_data(int reynolds) {
    GhiaData data;
    data.y_u = {1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0000};
    data.x_v = {1.0000,0.9688,0.9609,0.9531,0.9453,0.9063,0.8594,0.8047,0.5000,0.2344,0.2266,0.1563,0.0938,0.0781,0.0703,0.0625,0.0000};
    if (reynolds == 100) {
        data.u = {1.0000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717,0.0000};
        data.v = {0.0000,-0.05906,-0.07391,-0.08864,-0.10313,-0.16914,-0.22445,-0.24533,0.05454,0.17527,0.17507,0.16077,0.12317,0.10890,0.10091,0.09233,0.0000};
    } else if (reynolds == 400) {
        data.u = {1.0000,0.75837,0.68439,0.61756,0.55892,0.29093,0.16256,0.02135,-0.11477,-0.17119,-0.32726,-0.24299,-0.14612,-0.10338,-0.09266,-0.08186,0.0000};
        data.v = {0.0000,-0.12146,-0.15663,-0.19254,-0.22847,-0.23827,-0.44993,-0.38598,0.05186,0.30174,0.30203,0.28124,0.22965,0.20920,0.19713,0.18360,0.0000};
    } else if (reynolds == 1000) {
        data.u = {1.0000,0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,-0.06080,-0.10648,-0.27805,-0.38289,-0.29730,-0.22220,-0.20196,-0.18109,0.0000};
        data.v = {0.0000,-0.21388,-0.27669,-0.33714,-0.39188,-0.51550,-0.42665,-0.31966,0.02526,0.32235,0.33075,0.37095,0.32627,0.30353,0.29012,0.27485,0.0000};
    }
    return data;
}

static GhiaMetrics compare_with_ghia(const Result& result) {
    GhiaMetrics metrics;
    const GhiaData reference = ghia_data(result.reynolds);
    if (reference.u.empty()) {
        return metrics;
    }
    metrics.available = true;

    const int cells = result.cells;
    const int left = (cells - 1) / 2;
    const int right = cells / 2;
    std::vector<double> y_coordinates;
    std::vector<double> u_profile;
    std::vector<double> x_coordinates;
    std::vector<double> v_profile;
    y_coordinates.reserve(static_cast<std::size_t>(cells + 2));
    u_profile.reserve(static_cast<std::size_t>(cells + 2));
    x_coordinates.reserve(static_cast<std::size_t>(cells + 2));
    v_profile.reserve(static_cast<std::size_t>(cells + 2));

    y_coordinates.push_back(0.0);
    u_profile.push_back(0.0);
    for (int i = 0; i < cells; ++i) {
        y_coordinates.push_back(result.y[static_cast<std::size_t>(i)]);
        u_profile.push_back(0.5 * (result.u_center(i, left) + result.u_center(i, right)));
    }
    y_coordinates.push_back(1.0);
    u_profile.push_back(1.0);

    x_coordinates.push_back(0.0);
    v_profile.push_back(0.0);
    for (int j = 0; j < cells; ++j) {
        x_coordinates.push_back(result.x[static_cast<std::size_t>(j)]);
        v_profile.push_back(0.5 * (result.v_center(left, j) + result.v_center(right, j)));
    }
    x_coordinates.push_back(1.0);
    v_profile.push_back(0.0);

    double u_square_sum = 0.0;
    double v_square_sum = 0.0;
    for (std::size_t k = 0; k < reference.y_u.size(); ++k) {
        const double error = interpolate(y_coordinates, u_profile, reference.y_u[k]) - reference.u[k];
        u_square_sum += error * error;
        metrics.u_linf = std::isnan(metrics.u_linf) ? std::abs(error) : std::max(metrics.u_linf, std::abs(error));
    }
    for (std::size_t k = 0; k < reference.x_v.size(); ++k) {
        const double error = interpolate(x_coordinates, v_profile, reference.x_v[k]) - reference.v[k];
        v_square_sum += error * error;
        metrics.v_linf = std::isnan(metrics.v_linf) ? std::abs(error) : std::max(metrics.v_linf, std::abs(error));
    }
    metrics.u_l2 = std::sqrt(u_square_sum / static_cast<double>(reference.y_u.size()));
    metrics.v_l2 = std::sqrt(v_square_sum / static_cast<double>(reference.x_v.size()));

    if (result.reynolds == 100) {
        metrics.u_limit = 0.035;
        metrics.v_limit = 0.035;
    } else if (result.reynolds == 400) {
        metrics.u_limit = 0.10;
        metrics.v_limit = 0.13;
    } else {
        metrics.u_limit = 0.18;
        metrics.v_limit = 0.20;
    }
    metrics.passed = metrics.u_l2 <= metrics.u_limit && metrics.v_l2 <= metrics.v_limit;
    return metrics;
}

static std::string quality_label(const Result& result) {
    if (result.status == "converged" && result.ghia.available && result.ghia.passed) {
        return "converged_benchmark_pass";
    }
    if (result.status == "converged" && result.ghia.available) {
        return "converged_benchmark_needs_improvement";
    }
    if (result.status == "converged") {
        return "converged_no_benchmark";
    }
    if (result.ghia.available && result.ghia.passed) {
        return "benchmark_pass_not_converged";
    }
    return "needs_improvement";
}

static Result solve_case(
    int case_id,
    int cells,
    int reynolds,
    std::string scheme,
    std::string pressure_solver,
    const Config& config,
    const InitialState& initial = {}
) {
    scheme = lower(scheme);
    pressure_solver = upper(pressure_solver);
    const double h = config.length / static_cast<double>(cells);

    Field u(cells, cells + 1, 0.0);
    Field v(cells + 1, cells, 0.0);
    Field pressure(cells, cells, 0.0);
    if (initial.available
        && initial.u_face.rows == cells && initial.u_face.cols == cells + 1
        && initial.v_face.rows == cells + 1 && initial.v_face.cols == cells
        && initial.pressure.rows == cells && initial.pressure.cols == cells) {
        u = initial.u_face;
        v = initial.v_face;
        pressure = initial.pressure;
    }

    Result result;
    result.case_id = case_id;
    result.cells = cells;
    result.reynolds = reynolds;
    result.scheme = scheme;
    result.pressure_solver = pressure_solver;
    result.local_max_iterations = config.max_iterations;

    int consecutive_passes = 0;
    int consecutive_pressure_failures = 0;
    std::vector<double> stagnation_history;
    stagnation_history.reserve(static_cast<std::size_t>(config.stagnation_window));

    const auto start = std::chrono::steady_clock::now();
    for (int iteration = 1; iteration <= config.max_iterations; ++iteration) {
        const Field old_u = u;
        const Field old_v = v;

        auto [u_star, v_star, dt] = predict_velocity(
            u, v, pressure, reynolds, scheme, config
        );

        Field rhs = divergence(u_star, v_star, h, h, cells);
        for (double& value : rhs.values) {
            value /= dt;
        }

        auto [pressure_correction, poisson] = solve_pressure_poisson(
            rhs, h, pressure_solver, config
        );

        if (poisson.converged) {
            consecutive_pressure_failures = 0;
        } else {
            ++consecutive_pressure_failures;
            ++result.failed_pressure_solves;
        }

        u = u_star;
        v = v_star;
        for (int i = 0; i < cells; ++i) {
            for (int j = 1; j < cells; ++j) {
                u(i, j) -= dt * (pressure_correction(i, j) - pressure_correction(i, j - 1)) / h;
            }
        }
        for (int i = 1; i < cells; ++i) {
            for (int j = 0; j < cells; ++j) {
                v(i, j) -= dt * (pressure_correction(i, j) - pressure_correction(i - 1, j)) / h;
            }
        }
        for (std::size_t k = 0; k < pressure.values.size(); ++k) {
            pressure.values[k] += config.pressure_relaxation * pressure_correction.values[k];
        }
        remove_mean(pressure);

        const Residuals residuals = calculate_residuals(
            u, v, old_u, old_v, h, config
        );

        result.iterations = iteration;
        result.final_residuals = residuals;
        result.velocity_history.push_back(residuals.velocity_update_linf);
        result.divergence_linf_history.push_back(residuals.divergence_linf);
        result.divergence_l2_history.push_back(residuals.divergence_l2);
        result.mass_history.push_back(residuals.global_mass_imbalance);
        result.dt_history.push_back(dt);
        result.poisson_iteration_history.push_back(poisson.iterations);
        result.poisson_relative_history.push_back(poisson.relative_residual);
        result.poisson_converged_history.push_back(poisson.converged);

        if (!all_finite(u) || !all_finite(v) || !all_finite(pressure)
            || !std::isfinite(residuals.velocity_update_linf)
            || !std::isfinite(residuals.divergence_linf)
            || !std::isfinite(residuals.divergence_l2)) {
            result.status = "non_finite";
            break;
        }
        if (std::max({
                residuals.velocity_update_linf,
                residuals.divergence_linf,
                residuals.divergence_l2
            }) > config.diverged_limit) {
            result.status = "diverged";
            break;
        }
        if (consecutive_pressure_failures >= config.maximum_pressure_failures) {
            result.status = "pressure_not_converged";
            break;
        }

        const bool passed = poisson.converged
            && residuals.velocity_update_linf <= config.velocity_tolerance
            && residuals.divergence_linf <= config.divergence_linf_tolerance
            && residuals.divergence_l2 <= config.divergence_l2_tolerance
            && residuals.global_mass_imbalance <= 1e-12;

        if (iteration >= config.minimum_iterations && passed) {
            ++consecutive_passes;
        } else {
            consecutive_passes = 0;
        }
        result.consecutive_pass_count = consecutive_passes;
        if (consecutive_passes >= config.consecutive_passes) {
            result.status = "converged";
            break;
        }

        stagnation_history.push_back(residuals.velocity_update_linf);
        if (static_cast<int>(stagnation_history.size()) > config.stagnation_window) {
            stagnation_history.erase(stagnation_history.begin());
        }
        if (iteration > config.minimum_iterations
            && static_cast<int>(stagnation_history.size()) == config.stagnation_window) {
            const double first = stagnation_history.front();
            const double last = stagnation_history.back();
            const double relative_reduction = (first - last) / std::max(first, 1e-30);
            if (first > config.velocity_tolerance
                && relative_reduction < config.stagnation_minimum_reduction) {
                result.status = "stagnated";
                break;
            }
        }

        if (iteration % 1000 == 0) {
            std::cout << "      iter=" << iteration
                      << " vel=" << std::scientific << residuals.velocity_update_linf
                      << " div=" << residuals.divergence_linf
                      << " p=" << poisson.relative_residual << '\n';
        }
    }

    const auto end = std::chrono::steady_clock::now();
    result.runtime_seconds = std::chrono::duration<double>(end - start).count();
    result.u_face = std::move(u);
    result.v_face = std::move(v);
    result.pressure = std::move(pressure);
    build_cell_center_fields(result, config);
    result.ghia = compare_with_ghia(result);
    result.quality = quality_label(result);

    if (!result.poisson_iteration_history.empty()) {
        result.average_poisson_iterations = std::accumulate(
            result.poisson_iteration_history.begin(),
            result.poisson_iteration_history.end(),
            0.0
        ) / static_cast<double>(result.poisson_iteration_history.size());
        result.average_poisson_relative_residual = std::accumulate(
            result.poisson_relative_history.begin(),
            result.poisson_relative_history.end(),
            0.0
        ) / static_cast<double>(result.poisson_relative_history.size());
        int saturated = 0;
        for (const int pressure_iterations : result.poisson_iteration_history) {
            if (pressure_iterations >= config.poisson_max_iterations) {
                ++saturated;
            }
        }
        result.pressure_saturation_ratio = static_cast<double>(saturated)
                                         / static_cast<double>(result.poisson_iteration_history.size());
    }
    return result;
}

static std::string case_name(const Result& result) {
    std::ostringstream stream;
    stream << "case_" << std::setw(3) << std::setfill('0') << result.case_id
           << "_N" << result.cells
           << "_Re" << result.reynolds
           << '_' << result.scheme
           << '_' << result.pressure_solver
           << "_serial_cpp";
    return stream.str();
}

static void write_history(const Result& result, const Config& config) {
    const fs::path path = fs::path(config.data_directory) / (case_name(result) + "_history.csv");
    std::ofstream output(path);
    output << std::setprecision(12);
    output << "iter,velocity_update_linf,divergence_linf,divergence_l2,global_mass_imbalance,dt,poisson_iters,poisson_relative_residual,poisson_converged\n";
    for (std::size_t k = 0; k < result.velocity_history.size(); ++k) {
        output << (k + 1) << ','
               << result.velocity_history[k] << ','
               << result.divergence_linf_history[k] << ','
               << result.divergence_l2_history[k] << ','
               << result.mass_history[k] << ','
               << result.dt_history[k] << ','
               << result.poisson_iteration_history[k] << ','
               << result.poisson_relative_history[k] << ','
               << (result.poisson_converged_history[k] ? 1 : 0) << '\n';
    }
}

static void write_fields(const Result& result, const Config& config) {
    const fs::path path = fs::path(config.data_directory) / (case_name(result) + "_fields.csv");
    std::ofstream output(path);
    output << std::setprecision(12);
    output << "i,j,x,y,u,v,p,speed,vorticity\n";
    for (int i = 0; i < result.cells; ++i) {
        for (int j = 0; j < result.cells; ++j) {
            output << i << ',' << j << ','
                   << result.x[static_cast<std::size_t>(j)] << ','
                   << result.y[static_cast<std::size_t>(i)] << ','
                   << result.u_center(i, j) << ','
                   << result.v_center(i, j) << ','
                   << result.pressure(i, j) << ','
                   << result.speed(i, j) << ','
                   << result.vorticity(i, j) << '\n';
        }
    }
}

static void write_summary_header(std::ofstream& output) {
    output << "CaseID,Implementation,N,Re,Scheme,PressureSolver,Status,Quality,Iterations,LocalMaxIter,"
           << "FinalRu,FinalRv,FinalRcMass,FinalRcDiv,Runtime_s,AvgPoissonIterations,"
           << "AvgPoissonRelResidual,PressureSaturationRatio,HasGhia,ValidationPass,"
           << "Ghia_u_L2,Ghia_v_L2,Ghia_u_Linf,Ghia_v_Linf,Ghia_u_L2_Limit,Ghia_v_L2_Limit,"
           << "FinalVelocityLinf,FinalDivergenceL2,GlobalMassImbalance,FailedPressureSolves,ConsecutivePasses\n";
}

static void write_summary_row(std::ofstream& output, const Result& result) {
    output << std::setprecision(12)
           << result.case_id << ",serial_cpp,"
           << result.cells << ',' << result.reynolds << ','
           << result.scheme << ',' << result.pressure_solver << ','
           << result.status << ',' << result.quality << ','
           << result.iterations << ',' << result.local_max_iterations << ','
           << result.final_residuals.velocity_update_linf << ','
           << result.final_residuals.velocity_update_linf << ','
           << result.final_residuals.global_mass_imbalance << ','
           << result.final_residuals.divergence_linf << ','
           << result.runtime_seconds << ','
           << result.average_poisson_iterations << ','
           << result.average_poisson_relative_residual << ','
           << result.pressure_saturation_ratio << ','
           << (result.ghia.available ? 1 : 0) << ','
           << (result.ghia.passed ? 1 : 0) << ','
           << result.ghia.u_l2 << ',' << result.ghia.v_l2 << ','
           << result.ghia.u_linf << ',' << result.ghia.v_linf << ','
           << result.ghia.u_limit << ',' << result.ghia.v_limit << ','
           << result.final_residuals.velocity_update_linf << ','
           << result.final_residuals.divergence_l2 << ','
           << result.final_residuals.global_mass_imbalance << ','
           << result.failed_pressure_solves << ','
           << result.consecutive_pass_count << '\n';
}

static void configure_mode(
    const std::string& mode,
    Config& config,
    std::vector<int>& meshes,
    std::vector<int>& reynolds_numbers,
    std::vector<std::string>& schemes,
    std::vector<std::string>& pressure_solvers
) {
    const std::string normalized = lower(mode);
    if (normalized == "smoke") {
        meshes = {16};
        reynolds_numbers = {100};
        schemes = {"upwind"};
        pressure_solvers = {"RBGS"};
        config.max_iterations = 20;
        config.minimum_iterations = 100;
        config.poisson_max_iterations = 100;
        config.save_fields = false;
    } else if (normalized == "quick") {
        meshes = {24, 32};
        reynolds_numbers = {100};
        schemes = {"upwind", "central"};
        pressure_solvers = {"RBSOR"};
        config.velocity_tolerance = 1e-7;
    } else if (normalized == "medium") {
        meshes = {32};
        reynolds_numbers = {100, 400, 1000};
        schemes = {"upwind", "central"};
        pressure_solvers = {"RBSOR"};
        config.velocity_tolerance = 1e-7;
    } else if (normalized == "grid") {
        meshes = {16, 32, 64};
        reynolds_numbers = {100};
        schemes = {"central"};
        pressure_solvers = {"RBSOR"};
        config.velocity_tolerance = 1e-7;
    } else if (normalized == "full") {
        meshes = {32, 64, 128};
        reynolds_numbers = {100, 400, 1000};
        schemes = {"upwind", "central"};
        pressure_solvers = {"RBGS", "RBSOR"};
        config.velocity_tolerance = 1e-7;
    } else if (normalized != "single") {
        throw std::runtime_error("Unknown mode: " + mode + " (use smoke, single, quick, medium, grid, or full)");
    }
}

static void print_usage(const char* executable) {
    std::cout
        << "Usage:\n"
        << "  " << executable << " --mode smoke|single|quick|medium|grid|full\n"
        << "  " << executable << " --single --N 32 --Re 100 --scheme upwind --pressure RBSOR\n\n"
        << "Options:\n"
        << "  --no-fields                    Skip full field CSV output\n"
        << "  --strict                       Return a non-zero exit code when a case does not converge\n"
        << "  --maxIter VALUE                Maximum outer iterations\n"
        << "  --poisson-maxIter VALUE        Maximum pressure iterations\n"
        << "  --tol-velocity VALUE           Dimensionless velocity-update tolerance\n"
        << "  --tol-divergence VALUE         Dimensionless Linf divergence tolerance\n"
        << "  --tol-divergence-l2 VALUE      Dimensionless L2 divergence tolerance\n"
        << "  --poisson-tol VALUE            Relative pressure residual tolerance\n"
        << "  --alpha-u VALUE                Momentum relaxation factor\n"
        << "  --alpha-p VALUE                Pressure update relaxation factor\n"
        << "  --cfl VALUE                    Convective CFL limit\n"
        << "  --dt-max VALUE                 Maximum pseudo-time step\n"
        << "  --min-iterations VALUE         Minimum outer iterations\n"
        << "  --consecutive-passes VALUE     Required consecutive converged iterations\n";
}

int main(int argc, char** argv) {
    Config config;
    std::string mode = "quick";
    bool explicit_single = false;
    int single_cells = 32;
    int single_reynolds = 100;
    std::string single_scheme = "upwind";
    std::string single_pressure = "RBSOR";

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            auto require_value = [&]() -> std::string {
                if (i + 1 >= argc) {
                    throw std::runtime_error("Missing value for " + argument);
                }
                return argv[++i];
            };

            if (argument == "--help" || argument == "-h") {
                print_usage(argv[0]);
                return 0;
            } else if (argument == "--mode") {
                mode = require_value();
            } else if (argument == "--single") {
                mode = "single";
                explicit_single = true;
            } else if (argument == "--N") {
                single_cells = std::stoi(require_value());
                mode = "single";
                explicit_single = true;
            } else if (argument == "--Re") {
                single_reynolds = std::stoi(require_value());
                mode = "single";
                explicit_single = true;
            } else if (argument == "--scheme") {
                single_scheme = require_value();
                mode = "single";
                explicit_single = true;
            } else if (argument == "--pressure") {
                single_pressure = require_value();
                mode = "single";
                explicit_single = true;
            } else if (argument == "--implementation") {
                (void)require_value();
            } else if (argument == "--no-fields") {
                config.save_fields = false;
            } else if (argument == "--strict") {
                config.strict_exit = true;
            } else if (argument == "--maxIter") {
                config.max_iterations = std::stoi(require_value());
            } else if (argument == "--poisson-maxIter") {
                config.poisson_max_iterations = std::stoi(require_value());
            } else if (argument == "--tol-velocity") {
                config.velocity_tolerance = std::stod(require_value());
            } else if (argument == "--tol-divergence") {
                config.divergence_linf_tolerance = std::stod(require_value());
            } else if (argument == "--tol-divergence-l2") {
                config.divergence_l2_tolerance = std::stod(require_value());
            } else if (argument == "--poisson-tol") {
                config.poisson_relative_tolerance = std::stod(require_value());
            } else if (argument == "--alpha-u") {
                config.momentum_relaxation = std::stod(require_value());
            } else if (argument == "--alpha-p") {
                config.pressure_relaxation = std::stod(require_value());
            } else if (argument == "--cfl") {
                config.cfl = std::stod(require_value());
            } else if (argument == "--dt-max") {
                config.dt_max = std::stod(require_value());
            } else if (argument == "--min-iterations") {
                config.minimum_iterations = std::stoi(require_value());
            } else if (argument == "--consecutive-passes") {
                config.consecutive_passes = std::stoi(require_value());
            } else {
                throw std::runtime_error("Unknown argument: " + argument);
            }
        }

        std::vector<int> meshes;
        std::vector<int> reynolds_numbers;
        std::vector<std::string> schemes;
        std::vector<std::string> pressure_solvers;
        configure_mode(mode, config, meshes, reynolds_numbers, schemes, pressure_solvers);
        if (explicit_single || lower(mode) == "single") {
            meshes = {single_cells};
            reynolds_numbers = {single_reynolds};
            schemes = {lower(single_scheme)};
            pressure_solvers = {upper(single_pressure)};
        }

        fs::create_directories(config.data_directory);
        const fs::path summary_path = fs::path(config.data_directory) / ("study_summary_" + lower(mode) + ".csv");
        std::ofstream summary(summary_path);
        write_summary_header(summary);

        const int number_of_cases = static_cast<int>(
            meshes.size() * reynolds_numbers.size() * schemes.size() * pressure_solvers.size()
        );
        std::cout << "\nLID-DRIVEN CAVITY C++ PHASE 2 SOLVER\n"
                  << "Mode: " << mode << "\n"
                  << "Total simulations: " << number_of_cases << "\n"
                  << "Summary: " << summary_path << "\n\n";

        int case_id = 0;
        int failed_cases = 0;
        std::map<std::tuple<int, std::string, std::string>, InitialState> continuation_states;

        for (const int cells : meshes) {
            for (const std::string& scheme : schemes) {
                for (const std::string& pressure_solver : pressure_solvers) {
                    for (const int reynolds : reynolds_numbers) {
                        ++case_id;
                        const auto key = std::make_tuple(cells, lower(scheme), upper(pressure_solver));
                        InitialState initial;
                        const auto previous = continuation_states.find(key);
                        if (previous != continuation_states.end()) {
                            initial = previous->second;
                        }

                        std::cout << '[' << std::setw(3) << std::setfill('0') << case_id << std::setfill(' ')
                                  << "] N=" << cells
                                  << " Re=" << reynolds
                                  << " Scheme=" << scheme
                                  << " Pressure=" << pressure_solver << '\n';

                        Result result = solve_case(
                            case_id,
                            cells,
                            reynolds,
                            scheme,
                            pressure_solver,
                            config,
                            initial
                        );
                        write_summary_row(summary, result);
                        summary.flush();
                        write_history(result, config);
                        if (config.save_fields) {
                            write_fields(result, config);
                        }

                        if (result.status == "converged") {
                            continuation_states[key] = InitialState{
                                true,
                                result.u_face,
                                result.v_face,
                                result.pressure
                            };
                        } else {
                            ++failed_cases;
                        }

                        std::cout << "      status=" << result.status
                                  << " quality=" << result.quality
                                  << " iter=" << result.iterations << '/' << result.local_max_iterations
                                  << " vel=" << std::scientific << result.final_residuals.velocity_update_linf
                                  << " div=" << result.final_residuals.divergence_linf
                                  << " runtime=" << std::fixed << std::setprecision(2) << result.runtime_seconds << "s\n";
                        if (result.ghia.available) {
                            std::cout << "      Ghia L2: u=" << std::scientific << result.ghia.u_l2
                                      << " v=" << result.ghia.v_l2
                                      << " pass=" << (result.ghia.passed ? 1 : 0) << '\n';
                        }
                        std::cout << std::defaultfloat;
                    }
                }
            }
        }

        std::cout << "\nFinished. CSV outputs are in " << config.data_directory << '\n';
        if (config.strict_exit && failed_cases > 0) {
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        print_usage(argv[0]);
        return 1;
    }
}
