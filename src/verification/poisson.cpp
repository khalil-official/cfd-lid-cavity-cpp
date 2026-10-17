#include "lidcavity/poisson.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "lidcavity/operators.hpp"

namespace lidcavity {
namespace {

void validate_inputs(const Field2D& rhs, double dx, double dy, const PoissonOptions& options) {
    if (!(dx > 0.0) || !(dy > 0.0)) {
        throw std::invalid_argument("Grid spacing must be positive");
    }
    if (rhs.nx() < 3 || rhs.ny() < 3) {
        throw std::invalid_argument("Poisson grid requires at least 3x3 points");
    }
    if (options.max_iterations == 0 || options.check_every == 0) {
        throw std::invalid_argument("Iteration counts must be positive");
    }
    if (options.method == PoissonMethod::rbsor && !(options.omega > 0.0 && options.omega < 2.0)) {
        throw std::invalid_argument("SOR omega must be in (0, 2)");
    }
}

}  // namespace

double poisson_residual_linf(const Field2D& solution, const Field2D& rhs, double dx, double dy) {
    if (solution.nx() != rhs.nx() || solution.ny() != rhs.ny()) {
        throw std::invalid_argument("Solution and RHS shapes must match");
    }
    const Field2D laplacian = five_point_laplacian(solution, dx, dy);
    double residual = 0.0;
    for (std::size_t i = 1; i + 1 < rhs.ny(); ++i) {
        for (std::size_t j = 1; j + 1 < rhs.nx(); ++j) {
            residual = std::max(residual, std::abs(laplacian(i, j) - rhs(i, j)));
        }
    }
    return residual;
}

PoissonResult solve_poisson_dirichlet(
    const Field2D& rhs,
    double dx,
    double dy,
    const PoissonOptions& options
) {
    validate_inputs(rhs, dx, dy, options);

    Field2D solution(rhs.ny(), rhs.nx(), 0.0);
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double denominator = 2.0 * (dx2 + dy2);
    const double rhs_norm = std::max(linf_norm(rhs), 1.0);

    PoissonResult result{solution, 0, false, poisson_residual_linf(solution, rhs, dx, dy), 0.0};
    result.relative_residual = result.absolute_residual / rhs_norm;

    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        for (int color = 0; color < 2; ++color) {
            for (std::size_t i = 1; i + 1 < rhs.ny(); ++i) {
                const std::size_t first_j = 1U + ((i + 1U + static_cast<std::size_t>(color)) % 2U);
                for (std::size_t j = first_j; j + 1 < rhs.nx(); j += 2U) {
                    const double candidate =
                        ((solution(i + 1, j) + solution(i - 1, j)) * dx2
                         + (solution(i, j + 1) + solution(i, j - 1)) * dy2
                         - rhs(i, j) * dx2 * dy2) / denominator;
                    if (options.method == PoissonMethod::rbsor) {
                        solution(i, j) = (1.0 - options.omega) * solution(i, j)
                                       + options.omega * candidate;
                    } else {
                        solution(i, j) = candidate;
                    }
                }
            }
        }

        if (iteration == 1 || iteration % options.check_every == 0 || iteration == options.max_iterations) {
            result.absolute_residual = poisson_residual_linf(solution, rhs, dx, dy);
            result.relative_residual = result.absolute_residual / rhs_norm;
            if (result.absolute_residual <= options.absolute_tolerance
                || result.relative_residual <= options.relative_tolerance) {
                result.converged = true;
                result.iterations = iteration;
                result.solution = std::move(solution);
                return result;
            }
        }
    }

    result.iterations = options.max_iterations;
    result.solution = std::move(solution);
    return result;
}

std::string to_string(PoissonMethod method) {
    switch (method) {
        case PoissonMethod::rbgs: return "RBGS";
        case PoissonMethod::rbsor: return "RBSOR";
    }
    return "unknown";
}

}  // namespace lidcavity
