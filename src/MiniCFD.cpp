
#include "Grid.h"
#include "IO.h"
#include "Matrix.h"
#include "Solver.h"
#include "Vector.h"
#include <array>
#include <exception>
#include <iostream>
#include <memory>

using ScalarT = double;

/**
 * Helper classs to build blocked matrices using the existing SparseMatrix
 * class.
 */
template <typename T>
struct BlockMatBuilder {

    BlockMatBuilder<T>(size_t rows, size_t cols, size_t blockRows,
                       size_t blockCols)
        : rows(rows), cols(cols), blockRows(blockRows), blockCols(blockCols) {
        currentBlock.resize(blockRows * blockCols);
    }

    void finishBlock(size_t bi, size_t bj) {
        size_t idx;
        for (size_t i = 0; i < blockRows; i++) {
            for (size_t j = 0; j < blockCols; i++) {
                entries.push_back({bi * blockRows + i, bj * blockCols + j,
                                   currentBlock[idx]});
                idx++;
            }
        }
    }

    T& operator()(size_t i, size_t j) {
        assert("Access ouside of block bounds" && i >= 0 && i < blockRows &&
               j >= 0 && j < blockCols);
        return currentBlock[i * blockRows + j];
    }

    SparseMatrix<T> generate() const {
        return SparseMatrix<T>(rows, cols, entries);
    }

  private:
    size_t rows;
    size_t cols;
    size_t blockRows;
    size_t blockCols;
    std::vector<T> currentBlock;
    std::vector<SparseMatrixEntry<T>> entries;
};

template <typename T>
inline SparseMatrix<T> evalTransportJacobian(const VelocityField<T>& U) {
    size_t nPoints = U.getGrid()->getCellCount();
    BlockMatBuilder<T> builder(nPoints * 3, nPoints * 3, 3, 3);
    for (size_t i = 0; i < U.getWidth(); i++) {
        for (size_t j = 0; j < U.getHeight(); i++) {
            for (size_t k = 0; k < U.getDepth(); i++) {
                builder(0, 0) = U.dudx(i, j, k);
                builder(0, 1) = U.dudy(i, j, k);
                builder(0, 2) = U.dudz(i, j, k);
                builder(1, 0) = U.dvdx(i, j, k);
                builder(1, 1) = U.dvdy(i, j, k);
                builder(1, 2) = U.dvdz(i, j, k);
                builder(2, 0) = U.dwdx(i, j, k);
                builder(2, 1) = U.dwdy(i, j, k);
                builder(2, 2) = U.dwdz(i, j, k);
                builder.finishBlock(U.getGrid()->cellIndex(i, j, k),
                                    U.getGrid()->cellIndex(i, j, k));
            }
        }
    }
    return builder.generate();
}

template <typename T>
inline Vector<T> evalTransportEquation(const VelocityField<T>& U) {
    Vector<T> transportEq(U.getNumValues());
    size_t idx = 0;
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
                transportEq[idx++] = f_x;
                transportEq[idx++] = f_y;
                transportEq[idx++] = f_z;
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
    return startPos;
}

template <typename T>
inline VelocityField<T> semiLagrangianAdvection(const VelocityField<T>& U,
                                                double dt) {
    auto U_adv = U;
    auto cellSize = U.getCellSize();
    for (size_t i = 0; i < U.getWidth(); i++) {
        for (size_t j = 0; j < U.getHeight(); j++) {
            for (size_t k = 0; k < U.getDepth(); k++) {
                Vec3<T> pu{(i + 0.5) * cellSize, j * cellSize, k * cellSize};
                Vec3<T> pv{i * cellSize, (j + 0.5) * cellSize, k * cellSize};
                Vec3<T> pw{i * cellSize, j * cellSize, (k + 0.5) * cellSize};
                pu = traceBackward(pu, U, dt);
                pv = traceBackward(pv, U, dt);
                pw = traceBackward(pw, U, dt);
                U_adv.setLeftU(i, j, k, U.trilerp(pu).x);
                U_adv.setTopV(i, j, k, U.trilerp(pv).y);
                U_adv.setFrontW(i, j, k, U.trilerp(pw).z);
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
    auto size = U_adv.getGrid()->getCellCount();
    auto cellSize = U_adv.getGrid()->getCellSize();

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
    for (size_t k = 0; k < depth; k++) {
        for (size_t j = 0; j < height; j++) {
            for (size_t i = 0; i < width; i++) {
            
                // Right side
                auto U_div_ijk = U_adv.div(i, j, k);
                b[idx] = -(1.0 / dt) * U_div_ijk;

                // Number of non-BC neighboring cells
                int numRowEntries = 0;

                // Matrix coefficients
                if (i >= 1) {
                    assert(idx >= 1);
                    coeffs.push_back({idx, idx - 1, kNeighbor});
                    numRowEntries++;
                } else {
                    // Left BC
                }
                if (i < width - 1) {
                    coeffs.push_back({idx, idx + 1, kNeighbor});
                    numRowEntries++;
                } else {
                    // Right BC
                }
                if (j >= 1) {
                    assert(idx >= width);
                    coeffs.push_back({idx, idx - width, kNeighbor});
                    numRowEntries++;
                } else {
                    // Top BC
                }
                if (j < height - 1) {
                    coeffs.push_back({idx, idx + width, kNeighbor});
                    numRowEntries++;
                } else {
                    // Bottom BC
                }
                if (k >= 1) {
                    assert(idx >= width * height);
                    coeffs.push_back({idx, idx - width * height, kNeighbor});
                    numRowEntries++;
                } else {
                    // Front BC
                }
                if (k < depth - 1) {
                    coeffs.push_back({idx, idx + width * height, kNeighbor});
                    numRowEntries++;
                } else {
                    // Back BC
                }
                coeffs.push_back({idx, idx, numRowEntries * kMiddle});

                idx++;
            }
        }
    }
    SparseMatrix<T> A(size, size, coeffs);

    dumpMatrix(A);
    exit(0);

    std::cout << "Running PCG\n";
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

    for (size_t i = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++) {
            for (size_t k = 0; k < depth; k++) {
                auto u =
                    U_adv.getLeftU(i, j, k) -
                    dt * (p.getPressure(i, j, k) - p.getPressure(i - 1, j, k)) /
                        cellSize;
                auto v =
                    U_adv.getTopV(i, j, k) -
                    dt * (p.getPressure(i, j, k) - p.getPressure(i, j - 1, k)) /
                        cellSize;
                auto w =
                    U_adv.getFrontW(i, j, k) -
                    dt * (p.getPressure(i, j, k) - p.getPressure(i, j, k - 1)) /
                        cellSize;
                U_corr.setLeftU(i, j, k, u);
                U_corr.setTopV(i, j, k, v);
                U_corr.setFrontW(i, j, k, w);
            }
        }
    }
    return U_corr;
}

int main(int argc, char** argv) {
    // Command line parsing (later)

    auto N = 3;

    auto width = N;
    auto height = N;
    auto depth = N;

    // Set up grid
    auto grid = std::make_shared<Grid<ScalarT>>(width, height, depth, 1);

    // Create velovity and pressure fields
    auto U = std::make_unique<VelocityField<ScalarT>>(grid);
    auto p = std::make_unique<PressureField<ScalarT>>(grid);

    // Init values
    for (size_t i = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++) {
            for (size_t k = 0; k < depth; k++) {
                U->setLeftU(i, j, k, 0);
                U->setTopV(i, j, k, 0);
                U->setFrontW(i, j, k, 0);
                p->setPressure(i, j, k, 1);
                if (j == 1) {
                    U->setTopV(i, j, k, 5);
                }
            }
        }
    }

    double endTime = 10.0;
    double dt = 0.1;

    std::cout << "Initialization done!\n";
    write_to_file(*U, *p, "fields_0.txt");

    // Time loop
    unsigned step = 0;
    for (double t = 0; t < endTime; t += dt) {
        step++;
        // - Solve advection: implicit euler
        std::cout << "Solving advection equation\n";
        auto U_adv = solveAdvection(*U, dt);

        // - Pressure corection
        std::cout << "Solving pressure correction\n";
        auto p_new = solvePressureCorrection(U_adv, *p, dt);

        std::cout << "Applying pressure correction\n";
        auto U_corr = applyPressureCorrection(U_adv, p_new, dt);

        // - Write field
        write_to_file(U_corr, p_new, "fields_" + std::to_string(step) + ".txt");
        std::cout << "Time step complete. t = " << t << "\n";
    }

    // - Validation

    // - Visualization

    return 0;
}
