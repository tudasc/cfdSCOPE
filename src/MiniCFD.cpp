
#include "Grid.h"
#include "IO.h"
#include "Matrix.h"
#include "Solver.h"
#include "Util.h"
#include "Vector.h"

#include "spdlog/cfg/env.h"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"
#include <cxxopts.hpp>
#include <omp.h>

#include <array>
#include <exception>
#include <iostream>
#include <memory>

using ScalarT = double;

template <typename T>
inline Vector<T> evalTransportEquation(const VelocityField<T>& U) {
    Vector<T> transportEq(U.getNumValues());
    size_t idx = 0;
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < U.getWidth(); i++) {
        for (size_t j = 0; j < U.getHeight(); j++) {
            for (size_t k = 0; k < U.getDepth(); k++) {
                auto f_x = -U.getLeftU(i, j, k) * U.dudx(i, j, k) -
                           U.getTopV(i, j, k) * U.dudy(i, j, k) -
                           U.getFrontW(i, j, k) * U.dudz(i, j, k);
                auto f_y = -U.getLeftU(i, j, k) * U.dvdx(i, j, k) -
                           U.getTopV(i, j, k) * U.dvdy(i, j, k) -
                           U.getFrontW(i, j, k) * U.dvdz(i, j, k);
                auto f_z = -U.getLeftU(i, j, k) * U.dwdx(i, j, k) -
                           U.getTopV(i, j, k) * U.dwdy(i, j, k) -
                           U.getFrontW(i, j, k) * U.dwdz(i, j, k);
                transportEq[i * (U.getHeight()*U.getDepth()*3) + j*(U.getDepth()*3) + k*3 + 0] = f_x;
                transportEq[i * (U.getHeight()*U.getDepth()*3) + j*(U.getDepth()*3) + k*3 + 1] = f_y;
                transportEq[i * (U.getHeight()*U.getDepth()*3) + j*(U.getDepth()*3) + k*3 + 2] = f_z;
            }
        }
    }
    return transportEq;
}

template <typename T>
inline Vec3<T> traceBackward(Vec3<T> endPos, const VelocityField<T>& U,
                             float dt) {
    // Use 2nd order explicit Runge-Kutte scheme to trace velocity "particle"
    // back in time
    auto midPos = endPos - U.trilerp(endPos) * 0.5 * dt;
    auto startPos = endPos - U.trilerp(midPos) * dt;
    return U.getGrid()->getNearestInsidePos(startPos);
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
                    SPDLOG_TRACE("Cell at {}, {}, {}:  x comes from {}", pu.x,
                                 pu.y, pu.z, pu_origin.x);
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
inline PressureField<T> solvePressureCorrection(const VelocityField<T>& U_adv,
                                                const PressureField<T>& p,
                                                double dt) {
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

    //TODO Parallelize .push_back not suitable for openmp for
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

                SPDLOG_TRACE("Cell ({}, {}, {}): neighbors={}, div={}", i, j, k,
                             numRowEntries, U_div_ijk);
            }
        }
    }
    SparseMatrix<T> A(size, size, coeffs);

    SPDLOG_TRACE("Running PCG");
    Vector<T> p_new = pcg(A, b);
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
template<typename T>
inline VelocityField<T> applyForces(const VelocityField<T>& U,
                                                double dt) {

    auto width = U.getWidth();
    auto height = U.getHeight();
    auto depth = U.getDepth();

    auto cellSize = U.getCellSize();

    float windSpeed = 1;

    VelocityField<T> U_f = U;

    for (size_t i = 1; i < width; i++) { // Starting with i=1 to respect boundary condition
        for (size_t j = 0; j < height; j++) {
            for (size_t k = 0; k < depth; k++) {
                auto oldU = U_f.getLeftU(i, j, k);
                auto distToSurface = j * cellSize;
                auto maxForce = windSpeed - oldU;
                auto f = maxForce * std::exp(-2*distToSurface);           
                auto newU = oldU + dt * f;  // Assuming density rho = 1
                U_f.setLeftU(i, j, k, newU);                
            }
        }
    }
    return U_f;
}

int main(int argc, char** argv) {
    // Command line parsing
    cxxopts::Options options("MiniCFD",
                             "Simple CFD simulation for didactic purposes.");

    // clang-format off
    options.add_options()
        ("l,log-level", "Log level (trace, debug, info, warn, err, critical or off)", cxxopts::value<std::string>()->default_value("info"))
        ("d,domain-size", "Number of the simulation cells along all three axes", cxxopts::value<size_t>()->default_value("10"))
        ("c,cell-size", "Size of each simulation cell", cxxopts::value<ScalarT>()->default_value("1.0"))
        ("e,end-time", "Simulation duration (seconds)", cxxopts::value<double>()->default_value("1.0"))
        ("s,step-size", "Simulation step size (seconds)", cxxopts::value<double>()->default_value("0.05"))
        ("p,output-prefix", "Output file prefix", cxxopts::value<std::string>()->default_value("fields"))
        ("h,help", "Print usage")
    ;
    // clang-format on
    auto args = options.parse(argc, argv);
    if (args.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const std::string logLevel = args["log-level"].as<std::string>();
    if (logLevel == "trace") {
        spdlog::set_level(spdlog::level::trace);
    } else if (logLevel == "debug") {
        spdlog::set_level(spdlog::level::debug);
    } else if (logLevel == "info") {
        spdlog::set_level(spdlog::level::info);
    } else if (logLevel == "warn") {
        spdlog::set_level(spdlog::level::warn);
    } else if (logLevel == "err") {
        spdlog::set_level(spdlog::level::err);
    } else if (logLevel == "critical") {
        spdlog::set_level(spdlog::level::critical);
    } else if (logLevel == "off") {
        spdlog::set_level(spdlog::level::off);
    } else {
        spdlog::critical("Log level '{}' not recognized!");
        exit(-1);
    }

    // spdlog::cfg::load_env_levels();
    spdlog::info("Welcome to MiniCFD!");

    auto N = args["domain-size"].as<size_t>();
    auto width = N;
    auto height = N;
    auto depth = N;

    // Set up grid
    auto grid = std::make_shared<Grid<ScalarT>>(
        width, height, depth, args["cell-size"].as<ScalarT>());

    // Create velovity and pressure fields
    auto U = std::make_unique<VelocityField<ScalarT>>(grid);
    auto p = std::make_unique<PressureField<ScalarT>>(grid);

    // Init values
#pragma omp parallel
    {
#pragma omp master
        spdlog::info("Using {} OpenMP threads.", omp_get_num_threads());

#pragma omp for collapse(3)
        for (size_t i = 0; i < width; i++) {
            for (size_t j = 0; j < height; j++) {
                for (size_t k = 0; k < depth; k++) {
                    U->setLeftU(i, j, k, 0);
                    U->setTopV(i, j, k, 0);
                    U->setFrontW(i, j, k, 0);
                    p->setPressure(i, j, k, 1);
                }
            }
        }
    }

    double endTime = args["end-time"].as<double>();
    double dt = args["step-size"].as<double>();

    spdlog::info("Initialization done!");
    write_to_file(*U, *p, args["output-prefix"].as<std::string>() + "_0.txt");

    // Time loop
    unsigned step = 0;
    for (double t = 0; t < endTime; t += dt) {
        step++;

        SPDLOG_TRACE("Initial U_x field:\n{}",
                     dumpVectorComponent(U->getRawValues(), 0, 3));

        // - External forces
        spdlog::debug("Applying external forces");
        auto U_f = applyForces(*U, dt);
        SPDLOG_TRACE("U_x after forces:\n{}",
                     dumpVectorComponent(U_f.getRawValues(), 0, 3));

        // - Solve advection
        spdlog::debug("Solving advection equation");
        auto U_adv = solveAdvection(U_f, dt);

        SPDLOG_TRACE("U_x after advection:\n{}",
                     dumpVectorComponent(U_adv.getRawValues(), 0, 3));

        // - Pressure corection
        spdlog::debug("Solving pressure correction...");
        auto p_new = solvePressureCorrection(U_adv, *p, dt);

        SPDLOG_TRACE("Pressure field: {}\n", dumpVector(p_new.getRawValues()));

        spdlog::debug("Applying pressure correction...");
        auto U_corr = applyPressureCorrection(U_adv, p_new, dt);

        // - Write field
        write_to_file(U_corr, p_new,
                      args["output-prefix"].as<std::string>() + "_" +
                          std::to_string(step) + ".txt");
        spdlog::info("Time step complete. t = {:.5f}", t);

        *U = U_corr;
        *p = p_new;
    }

    return 0;
}
