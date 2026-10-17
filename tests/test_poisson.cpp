#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "lidcavity/poisson.hpp"

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

struct ErrorResult {
    double l2 = 0.0;
    lidcavity::PoissonResult solve;
};

ErrorResult run_case(std::size_t n, lidcavity::PoissonMethod method) {
    using namespace lidcavity;
    const double h = 1.0 / static_cast<double>(n - 1);
    Field2D rhs(n, n, 0.0);
    Field2D exact(n, n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const double x = static_cast<double>(j) * h;
            const double y = static_cast<double>(i) * h;
            exact(i, j) = std::sin(pi * x) * std::sin(pi * y);
            rhs(i, j) = -2.0 * pi * pi * exact(i, j);
        }
    }

    PoissonOptions options;
    options.method = method;
    options.max_iterations = 30000;
    options.absolute_tolerance = 1e-11;
    options.relative_tolerance = 1e-11;
    options.omega = 1.7;
    options.check_every = 10;

    PoissonResult solve = solve_poisson_dirichlet(rhs, h, h, options);
    require(solve.converged, "manufactured Poisson case converged");

    long double sum = 0.0L;
    std::size_t count = 0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        for (std::size_t j = 1; j + 1 < n; ++j) {
            const long double error = static_cast<long double>(solve.solution(i, j) - exact(i, j));
            sum += error * error;
            ++count;
        }
    }

    return {std::sqrt(static_cast<double>(sum / static_cast<long double>(count))), std::move(solve)};
}

}  // namespace

int main() {
    using namespace lidcavity;
    // A quadratic manufactured solution is exact under the five-point stencil.
    // Rectangular grids exercise both starting parities and unequal spacings.
    for (const auto shape : {std::pair<std::size_t, std::size_t>{8, 11}, {9, 10}}) {
        const auto ny = shape.first;
        const auto nx = shape.second;
        const double dx = 2.0 / static_cast<double>(nx - 1);
        const double dy = 1.0 / static_cast<double>(ny - 1);
        Field2D rhs(ny, nx, 0.0);
        for (std::size_t i = 0; i < ny; ++i) {
            for (std::size_t j = 0; j < nx; ++j) {
                const double x = static_cast<double>(j) * dx;
                const double y = static_cast<double>(i) * dy;
                rhs(i, j) = -2.0 * (y * (1.0 - y) + x * (2.0 - x));
            }
        }
        for (const auto method : {PoissonMethod::rbgs, PoissonMethod::rbsor}) {
            PoissonOptions options;
            options.method = method;
            options.absolute_tolerance = 1e-11;
            options.relative_tolerance = 1e-11;
            const auto result = solve_poisson_dirichlet(rhs, dx, dy, options);
            require(result.converged, "rectangular Poisson case converged");
            for (std::size_t i = 0; i < ny; ++i) {
                for (std::size_t j = 0; j < nx; ++j) {
                    const double x = static_cast<double>(j) * dx;
                    const double y = static_cast<double>(i) * dy;
                    require(std::abs(result.solution(i, j) - x * (2.0 - x) * y * (1.0 - y)) < 1e-9,
                            "rectangular solution matches exact quadratic including boundaries");
                }
            }
        }
    }
    const std::vector<std::size_t> grids{17, 33, 65};
    std::vector<double> errors;

    for (const std::size_t n : grids) {
        const ErrorResult result = run_case(n, PoissonMethod::rbsor);
        errors.push_back(result.l2);
        std::cout << "N=" << n << " L2=" << result.l2
                  << " iterations=" << result.solve.iterations
                  << " relres=" << result.solve.relative_residual << '\n';
    }

    const double order_1 = std::log(errors[0] / errors[1]) / std::log(2.0);
    const double order_2 = std::log(errors[1] / errors[2]) / std::log(2.0);
    require(order_1 > 1.8 && order_2 > 1.8, "Poisson discretization shows approximately second-order convergence");

    const ErrorResult rbgs = run_case(33, PoissonMethod::rbgs);
    const ErrorResult rbsor = run_case(33, PoissonMethod::rbsor);
    double method_difference = 0.0;
    for (std::size_t i = 0; i < rbgs.solve.solution.ny(); ++i) {
        for (std::size_t j = 0; j < rbgs.solve.solution.nx(); ++j) {
            method_difference = std::max(
                method_difference,
                std::abs(rbgs.solve.solution(i, j) - rbsor.solve.solution(i, j))
            );
        }
    }
    require(method_difference < 1e-8, "RBGS and RBSOR converge to the same discrete solution");

    std::cout << "Observed orders: " << order_1 << ", " << order_2
              << "; RBGS/RBSOR difference=" << method_difference << '\n';
    return EXIT_SUCCESS;
}
