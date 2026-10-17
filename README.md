# Lid-Driven Cavity Flow Solver in C++

Development edition based on [Ahmed Kandil's original project](https://github.com/Kandil2001/LidCavity_CPP).
Original authorship, Git history, and the MIT license are retained.

## Local improvement and Windows build

Red-black Poisson sweeps visit only the selected color with a stride of two,
removing the skipped half of the inner-loop iterations while preserving update
order. This is applied to both the production solver and the verification solver.
CMake selects MSVC warning flags on Windows; passing GCC's `-Wextra` to MSVC
previously prevented compilation.

From a PowerShell session with CMake available:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

Verified on 2026-10-18 using MSVC 19.44 and CMake 3.31.6. All four CTest tests
pass, including new rectangular manufactured Poisson cases with both grid
parities and both RBGS/RBSOR. The canonical N=32, Re=100 regression retained
2,086 iterations and identical CSV numerical metrics. One local Release run
decreased from 2.80 s to 1.65 s; this single-run timing is indicative only.

---

<p align="center">
  <img src="https://img.shields.io/badge/Status-Phase%202%20verified-brightgreen.svg" alt="Phase 2 verified">
  <img src="https://img.shields.io/badge/C%2B%2B-17-blue.svg" alt="C++17">
  <img src="https://img.shields.io/badge/Grid-staggered%20MAC-blueviolet.svg" alt="Staggered MAC grid">
  <a href="https://github.com/Kandil2001/LidCavity_CPP/actions/workflows/ci.yml">
    <img src="https://github.com/Kandil2001/LidCavity_CPP/actions/workflows/ci.yml/badge.svg" alt="Build and smoke test">
  </a>
  <a href="https://github.com/Kandil2001/LidCavity_CPP/actions/workflows/verification.yml">
    <img src="https://github.com/Kandil2001/LidCavity_CPP/actions/workflows/verification.yml/badge.svg" alt="Numerical verification">
  </a>
  <img src="https://img.shields.io/badge/License-MIT-lightgrey.svg" alt="MIT License">
</p>

A serial C++17 solver for the two-dimensional incompressible lid-driven cavity benchmark.

The production solver uses a staggered Marker-and-Cell arrangement, a pseudo-transient projection method, compatible pressure-gradient and divergence operators, and convergence-aware termination. The repository also contains manufactured Poisson verification, discrete-operator tests, sanitizer builds, Ghia centerline comparisons, and configurable parameter studies.

## What is verified

The current Phase 2 regression set has been exercised with:

- the canonical `N=32`, `Re=100`, upwind, RBSOR case;
- all six `N=32` combinations at `Re=100`, `400`, and `1000` with upwind and central convection using RBSOR;
- the `N=16`, `32`, and `64`, `Re=100`, central, RBSOR grid sequence;
- RBGS/RBSOR manufactured Poisson tests;
- the compatibility test between the selected divergence, gradient, and Laplacian operators;
- GCC and Clang builds in Debug and Release configurations;
- AddressSanitizer and UndefinedBehaviorSanitizer checks.

The canonical strict case converges automatically rather than stopping at a configured iteration limit.

## Numerical method

The solver advances the nondimensional incompressible Navier–Stokes equations using:

1. a staggered MAC grid;
2. an explicit pseudo-time momentum predictor;
3. upwind or central convection differencing;
4. a pressure-correction Poisson equation;
5. RBGS or RBSOR pressure iteration;
6. velocity correction using pressure gradients located consistently with the staggered velocity components;
7. convergence checks based on velocity updates, local divergence, global mass balance, and pressure convergence.

Velocity components are stored on cell faces and pressure is stored at cell centers. Cell-centered velocity values are reconstructed for CSV output and post-processing.

## Convergence contract

A production case is reported as `converged` only when all of the following remain satisfied for the configured number of consecutive iterations:

- dimensionless velocity-update `Linf` residual;
- dimensionless divergence `Linf` residual;
- dimensionless divergence `L2` residual;
- global boundary mass imbalance;
- successful pressure-Poisson convergence.

Possible terminal states include:

- `converged`
- `max_iterations`
- `pressure_not_converged`
- `stagnated`
- `diverged`
- `non_finite`

Use `--strict` when a script or CI job must fail if any requested case does not converge.

## Quick start

On Linux, WSL, or a Linux HPC node:

```bash
bash scripts/run_single.sh
```

This builds the solver and runs the canonical strict case:

```text
N = 32
Re = 100
scheme = upwind
pressure solver = RBSOR
```

A successful run ends with `status=converged` and writes CSV files to `results/data/`.

## Available runs

```bash
bash scripts/run_smoke_test.sh   # compilation and output check only
bash scripts/run_single.sh       # canonical converged regression
bash scripts/run_quick.sh        # four fast Re=100 cases
bash scripts/run_medium.sh       # six N=32 cases at Re=100/400/1000
bash scripts/run_grid.sh         # N=16/32/64 grid sequence at Re=100
bash scripts/run_re1000.sh       # converged N=32, Re=1000 central case
bash scripts/run_full.sh         # complete 36-case configuration
```

The full mode includes RBGS and grids up to `N=128`; it is intended for a workstation or HPC node.

## Direct command-line use

```bash
bash scripts/build.sh

bin/lid_cavity \
  --single \
  --N 32 \
  --Re 100 \
  --scheme upwind \
  --pressure RBSOR \
  --strict
```

Important numerical options include:

```text
--maxIter
--poisson-maxIter
--tol-velocity
--tol-divergence
--tol-divergence-l2
--poisson-tol
--alpha-u
--alpha-p
--cfl
--dt-max
--min-iterations
--consecutive-passes
```

Run `bin/lid_cavity --help` for the complete interface.

## CMake and tests

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The CTest suite includes:

- discrete operator verification;
- manufactured Poisson verification and grid refinement;
- convergence-state logic;
- the canonical production-solver regression.

For sanitizer checks:

```bash
cmake -S . -B build/sanitized \
  -DCMAKE_BUILD_TYPE=Debug \
  -DLIDCAVITY_ENABLE_SANITIZERS=ON
cmake --build build/sanitized --parallel
ctest --test-dir build/sanitized --output-on-failure
```

## Output files

Each case writes:

```text
results/data/study_summary_<mode>.csv
results/data/case_<id>_..._history.csv
results/data/case_<id>_..._fields.csv
```

The summary separates:

- execution status;
- iterative convergence;
- pressure convergence;
- divergence metrics;
- runtime;
- Ghia benchmark agreement.

The history file stores iteration-by-iteration convergence information. The field file contains cell-center coordinates, velocity, pressure, speed, and vorticity.

## Ghia benchmark comparison

For `Re=100`, `400`, and `1000`, the solver compares:

- horizontal velocity `u(y)` on the vertical centerline;
- vertical velocity `v(x)` on the horizontal centerline.

The code reports `L2` and `Linf` errors. These are reference-benchmark comparisons, not experimental validation.

## Repository structure

```text
src/lid_cavity.cpp              production staggered-grid solver
include/lidcavity/              reusable verification interfaces
src/verification/               operator, Poisson, and convergence components
tests/                          CTest verification programs
scripts/                        build and run helpers
postprocess/                    Python plotting scripts
docs/                           method and verification notes
results/data/                   generated CSV output
.github/workflows/              smoke and numerical-verification CI
```

## Scope and limitations

This is an educational CFD and numerical-verification project, not a production CFD package.

Current limitations include:

- serial CPU execution only;
- explicit pseudo-time momentum advancement;
- no multigrid pressure solver;
- no formal experimental validation;
- the exhaustive 36-case study can be computationally expensive, especially with RBGS and `N=128`.

Parallel C++, MPI, OpenMP, CUDA, MATLAB, and Python comparisons belong to the separate solver-comparison project and are intentionally not developed in this repository.

## Reference

Ghia, U., Ghia, K. N., & Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method*. Journal of Computational Physics, 48(3), 387–411.

## Author

Ahmed Kandil — [Portfolio](https://kandil2001.github.io/) · [LinkedIn](https://www.linkedin.com/in/ahmed-kandil03/) · [ORCID](https://orcid.org/0009-0007-2724-4565)

Released under the [MIT License](LICENSE).
