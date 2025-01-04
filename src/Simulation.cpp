#include "Simulation.h"

#include "Grid.h"
#include "IO.h"
#include "Matrix.h"
#include "Solver.h"
#include "Vector.h"

#include <omp.h>

template <typename T>
inline Vec3<T> traceBackward(Vec3<T> endPos, const VelocityField<T>& U,
                             float dt) {
    // Use 4th order explicit Runge-Kutte scheme to trace velocity "particle"
    // back in time
    auto k1 = U.trilerp(endPos);
    auto k2 = U.trilerp(endPos - k1 * (dt * 0.5));
    auto k3 = U.trilerp(endPos - k2 * (dt * 0.5));
    auto k4 = U.trilerp(endPos - k3 * dt);

    auto res = endPos - (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6.0);
    return U.getGrid()->getNearestInsidePos(res);
}

template <typename T>
inline VelocityField<T> semiLagrangianAdvection(const VelocityField<T>& U,
                                                double dt) {
    auto U_adv = U;
    auto cellSize = U.getCellSize();
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < U.getWidth(); i++) {
        for (size_t j = 0; j < U.getHeight(); j++) {
            for (size_t k = 0; k < U.getDepth(); k++) {
                // u is on boundary for i = 0
                if (i > 0) {
                    Vec3<T> pu{i * cellSize, (j + 0.5) * cellSize,
                               (k + 0.5) * cellSize};
                    auto pu_origin = traceBackward(pu, U, dt);
                    U_adv.setLeftU(i, j, k, U.trilerp(pu_origin).x);
                    // SPDLOG_TRACE("Cell at {}, {}, {}:  x comes from {}",
                    // pu.x,
                    //              pu.y, pu.z, pu_origin.x);
                }
                // v is on boundary for j = 0
                if (j > 0) {
                    Vec3<T> pv{(i + 0.5) * cellSize, j * cellSize,
                               (k + 0.5) * cellSize};
                    auto pv_origin = traceBackward(pv, U, dt);
                    U_adv.setTopV(i, j, k, U.trilerp(pv_origin).y);
                }
                // w is on boundary for k = 0
                if (k > 0) {
                    Vec3<T> pw{(i + 0.5) * cellSize, (j + 0.5) * cellSize,
                               k * cellSize};
                    auto pw_origin = traceBackward(pw, U, dt);
                    U_adv.setFrontW(i, j, k, U.trilerp(pw_origin).z);
                }
            }
        }
    }
    return U_adv;
}

template <typename T>
inline VelocityField<T> solveAdvection(const VelocityField<T>& U0, double dt) {
    auto U_adv = semiLagrangianAdvection(U0, dt);
    return U_adv;
}

template <typename T>
inline PressureField<T> solvePressureCorrection(
    const VelocityField<T>& U_adv, const PressureField<T>& p, double dt,
    PreconditionerType preconditioner = PreconditionerType::DIC) {
    auto& grid = *U_adv.getGrid();

    auto size = grid.getCellCount();
    auto cellSize = grid.getCellSize();

    auto width = U_adv.getWidth();
    auto height = U_adv.getHeight();
    auto depth = U_adv.getDepth();

    Vector<T> b(size);
    Vector<T> p_laplace(size);

    // Scaling factors for laplace stencil operator
    float kNeighbor = -1.0 / (cellSize * cellSize);
    float kMiddle = 1.0 / (cellSize * cellSize);

    // Note: Currently assuming Neumann for all pressure BCs
    // TODO: Implement Dirichlet BCs

    size_t idx = 0;
    std::vector<SparseMatrixEntry<T>> coeffs;

    // TODO Parallelize .push_back not suitable for openmp for
    for (size_t k = 0; k < depth; k++) {
        for (size_t j = 0; j < height; j++) {
            for (size_t i = 0; i < width; i++) {

                size_t idx = grid.cellIndex(i, j, k);

                // Right side
                auto U_div_ijk = U_adv.div(i, j, k);
                b[idx] = -(1.0 / dt) * U_div_ijk;

                // Number of non-BC neighboring cells
                int numRowEntries = 0;

                // Matrix coefficients
                if (i >= 1) {
                    assert(idx >= 1);
                    coeffs.push_back(
                        {idx, grid.cellIndex(i - 1, j, k), kNeighbor});
                    numRowEntries++;
                } else {
                    // Left BC
                }
                if (i < width - 1) {
                    coeffs.push_back(
                        {idx, grid.cellIndex(i + 1, j, k), kNeighbor});
                    numRowEntries++;
                } else {
                    // Right BC
                }
                if (j >= 1) {
                    assert(idx >= width);
                    coeffs.push_back(
                        {idx, grid.cellIndex(i, j - 1, k), kNeighbor});
                    numRowEntries++;
                } else {
                    // Top BC
                }
                if (j < height - 1) {
                    coeffs.push_back(
                        {idx, grid.cellIndex(i, j + 1, k), kNeighbor});
                    numRowEntries++;
                } else {
                    // Bottom BC
                }
                if (k >= 1) {
                    assert(idx >= width * height);
                    coeffs.push_back(
                        {idx, grid.cellIndex(i, j, k - 1), kNeighbor});
                    numRowEntries++;
                } else {
                    // Front BC
                }
                if (k < depth - 1) {
                    coeffs.push_back(
                        {idx, grid.cellIndex(i, j, k + 1), kNeighbor});
                    numRowEntries++;
                } else {
                    // Back BC
                }
                coeffs.push_back({idx, idx, numRowEntries * kMiddle});

                // SPDLOG_TRACE("Cell ({}, {}, {}): neighbors={}, div={}", i, j,
                // k,
                //              numRowEntries, U_div_ijk);
            }
        }
    }
    SparseMatrix<T> A(size, size, coeffs);

    SPDLOG_TRACE("Running PCG");
    std::unique_ptr<Preconditioner<T>> pre;
    switch (preconditioner) {
    case PreconditionerType::NONE:
        pre = std::make_unique<Preconditioner<T>>();
        break;
    case PreconditionerType::JACOBI:
        pre = std::make_unique<JacobiPreconditioner<T>>(A);
        break;
    case PreconditionerType::DIC:
        pre = std::make_unique<DICPreconditioner<T>>(A);
        break;
    default:
        assert("Unknown preconditioner type" && false);
    }

    Vector<T> p_new = pcg(A, b, *pre, p.getRawValues());
    return {p.getGrid(), p_new};
}

template <typename T>
inline VelocityField<T> applyPressureCorrection(const VelocityField<T>& U_adv,
                                                const PressureField<T>& p,
                                                double dt) {

    auto width = U_adv.getWidth();
    auto height = U_adv.getHeight();
    auto depth = U_adv.getDepth();

    auto cellSize = p.getCellSize();

    VelocityField<T> U_corr = U_adv;

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++) {
            for (size_t k = 0; k < depth; k++) {
                // u is on boundary for i = 0
                if (i > 0) {
                    auto u = U_adv.getLeftU(i, j, k) -
                             dt *
                                 (p.getPressure(i, j, k) -
                                  p.getPressure(i - 1, j, k)) /
                                 cellSize;
                    U_corr.setLeftU(i, j, k, u);
                }
                // v is on boundary for j = 0
                if (j > 0) {
                    auto v = U_adv.getTopV(i, j, k) -
                             dt *
                                 (p.getPressure(i, j, k) -
                                  p.getPressure(i, j - 1, k)) /
                                 cellSize;
                    U_corr.setTopV(i, j, k, v);
                }
                // w is on boundary for k = 0
                if (k > 0) {
                    auto w = U_adv.getFrontW(i, j, k) -
                             dt *
                                 (p.getPressure(i, j, k) -
                                  p.getPressure(i, j, k - 1)) /
                                 cellSize;
                    U_corr.setFrontW(i, j, k, w);
                }
            }
        }
    }
    return U_corr;
}

/**
 * Applying a constant force (e.g. wind) on the fluid surface.
 */
template <typename T>
inline VelocityField<T> applyForces(const VelocityField<T>& U, double dt) {

    auto width = U.getWidth();
    auto height = U.getHeight();
    auto depth = U.getDepth();

    auto cellSize = U.getCellSize();

    float windSpeed = 10 * cellSize;

    VelocityField<T> U_f = U;

    for (size_t i = 1; i < width;
         i++) { // Starting with i=1 to respect boundary condition
        for (size_t j = 0; j < height; j++) {
            for (size_t k = 0; k < depth; k++) {
                auto oldU = U_f.getLeftU(i, j, k);
                auto distToSurface = j * cellSize;
                auto maxForce = windSpeed - oldU;
                auto f = maxForce * std::exp(-2 * distToSurface);
                auto newU = oldU + dt * f; // Assuming density rho = 1
                U_f.setLeftU(i, j, k, newU);
            }
        }
    }
    return U_f;
}

SimulationOutput simulate(const SimulationConfig& cfg) {
    // Set up grid
    auto grid = std::make_shared<Grid<ScalarT>>(cfg.width, cfg.height,
                                                cfg.depth, cfg.cellSize);

    // Create velovity and pressure fields
    auto U = std::make_unique<VelocityField<ScalarT>>(grid);
    auto p = std::make_unique<PressureField<ScalarT>>(grid);

    // Init values
#pragma omp parallel
    {
#pragma omp master
        spdlog::info("Using {} OpenMP threads.", omp_get_num_threads());

#pragma omp for collapse(3)
        for (size_t i = 0; i < cfg.width; i++) {
            for (size_t j = 0; j < cfg.height; j++) {
                for (size_t k = 0; k < cfg.depth; k++) {
                    U->setLeftU(i, j, k, 0);
                    U->setTopV(i, j, k, 0);
                    U->setFrontW(i, j, k, 0);
                    p->setPressure(i, j, k, 1);
                }
            }
        }
    }

    spdlog::info("Initialization done!");
    if (!cfg.disableFileOutput) {
        write_to_file(*U, *p, cfg.outputPrefix + "_0.txt");
    }

    VelocityField<ScalarT> U_prev = *U;

    // Time loop
    unsigned step = 0;
    double dt = cfg.stepSize;
    for (double t = 0; t < cfg.endTime; t += dt) {
        step++;

        // SPDLOG_TRACE("Initial U_x field:\n{}",
        //              dumpVectorComponent(U->getRawValues(), 0, 3));

        // - External forces
        spdlog::debug("Applying external forces");
        auto U_f = applyForces(*U, dt);
        // SPDLOG_TRACE("U_x after forces:\n{}",
        //              dumpVectorComponent(U_f.getRawValues(), 0, 3));

        // - Solve advection
        spdlog::debug("Solving advection equation");
        auto U_adv = solveAdvection(U_f, dt);

        // SPDLOG_TRACE("U_x after advection:\n{}",
        //              dumpVectorComponent(U_adv.getRawValues(), 0, 3));

        // - Pressure corection
        spdlog::debug("Solving pressure correction...");
        auto p_new =
            solvePressureCorrection(U_adv, *p, dt, cfg.preconditionerType);

        // SPDLOG_TRACE("Pressure field: {}\n",
        // dumpVector(p_new.getRawValues()));

        spdlog::debug("Applying pressure correction...");
        auto U_corr = applyPressureCorrection(U_adv, p_new, dt);

        // - Write field
        if (!cfg.disableFileOutput) {
            write_to_file(U_corr, p_new,
                          cfg.outputPrefix + "_" + std::to_string(step) +
                              ".txt");
        }
        spdlog::info("Time step complete. t = {:.5f}", t);

        *U = U_corr;
        *p = p_new;

        auto dif = norm(U->getRawValues() - U_prev.getRawValues());
        spdlog::debug("Difference in U: {:.5f}", dif);
        U_prev = *U;
    }

    return SimulationOutput{std::move(U), std::move(p)};
}
