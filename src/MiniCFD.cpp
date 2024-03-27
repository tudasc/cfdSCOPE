
#include "IO.h"
#include "Grid.h"
#include "Matrix.h"
#include "Solver.h"
#include <array>
#include <iostream>
#include <exception>
#include <memory>

using ScalarT = double;


/**
* Helper classs to build blocked matrices using the existing SparseMatrix class.
*/
template<typename T> 
struct BlockMatBuilder {

    BlockMatBuilder<T>(size_t rows, size_t cols, size_t blockRows, size_t blockCols) :  rows(rows), cols(cols), blockRows(blockRows), blockCols(blockCols) {
        currentBlock.resize(blockRows * blockCols);
    }

    void finishBlock(size_t bi, size_t bj) { 
        size_t idx;
        for (size_t i = 0; i < blockRows; i++) {
            for (size_t j = 0; j < blockCols; i++) {
                entries.push_back({bi * blockRows + i, bj * blockCols + j, currentBlock[idx]});
                idx++;
            }
        }
    }

    T& operator()(size_t i, size_t j) {
        assert("Access ouside of block bounds" && i >= 0 && i < blockRows && j >= 0 && j < blockCols);
        return currentBlock[i*blockRows+j]; 
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

template<typename T>
inline SparseMatrix<T> evalTransportJacobian(const VelocityField<T>& v) {
    size_t nPoints = v.getGrid()->getCellCount();
    BlockMatBuilder<T> builder(nPoints * 3, nPoints * 3, 3, 3);
    for (size_t i = 0; i < v.getWidth(); i++) {
        for (size_t j = 0; j < v.getHeight(); i++) {
            for (size_t k = 0; k < v.getDepth(); i++) {
                builder(0, 0) = v.dudx(i, j, k);
                builder(0, 1) = v.dudy(i, j, k);
                builder(0, 2) = v.dudz(i, j, k);
                builder(1, 0) = v.dvdx(i, j, k);
                builder(1, 1) = v.dvdy(i, j, k);
                builder(1, 2) = v.dvdz(i, j, k);
                builder(2, 0) = v.dwdx(i, j, k);
                builder(2, 1) = v.dwdy(i, j, k);
                builder(2, 2) = v.dwdz(i, j, k);
                builder.finishBlock(v.getGrid()->cellIndex(i, j, k), v.getGrid()->cellIndex(i, j, k));
            }
        }
    }
    return builder.generate();
}

template<typename T>
inline Vector<T> evalTransportEquation(const VelocityField<T>& v) {
    Vector<T> transportEq(v.getNumValues());
    size_t idx = 0;
    for (size_t i = 0; i < v.getWidth(); i++) {
        for (size_t j = 0; j < v.getHeight(); i++) {
            for (size_t k = 0; k < v.getDepth(); i++) {
                auto f_x = - v.getLeftU(i, j, k) * v.dudx(i, j, k) - v.getTopV(i, j, k) * v.dudy(i, j, k) - v.getFrontW(i, j, k) * v.dudz(i,j,k);
                auto f_y = - v.getLeftU(i, j, k) * v.dvdx(i, j, k) - v.getTopV(i, j, k) * v.dvdy(i, j, k) - v.getFrontW(i, j, k) * v.dvdz(i,j,k);
                auto f_z = - v.getLeftU(i, j, k) * v.dwdx(i, j, k) - v.getTopV(i, j, k) * v.dwdy(i, j, k) - v.getFrontW(i, j, k) * v.dwdz(i,j,k);
                transportEq[idx++] = f_x;
                transportEq[idx++] = f_y;
                transportEq[idx++] = f_z;
            }
        }
    }
    return transportEq;
}

template<typename T>
inline VelocityField<T> solveAdvection(const VelocityField<T>& v0, double dt) {
    auto J = evalTransportJacobian(v0);

    // Assemble full matrix for implicit Euler: I - dt*J
    auto A = J * (-dt);
    for (size_t ij = 0; ij < A.getCols(); ij++) {
        A(ij,ij) = A(ij,ij) + 1;
    }
    
    auto b = evalTransportEquation(v0) * dt;

    auto v_adv = pcg(A, b);

    return {v0.getGrid(), v_adv};
}

template<typename T>
inline void applyPressureCorrection(const VelocityField<T>& v_adv, const PressureField<T>& p) {
    auto size = v_adv.getGrid()->getCellCount();
    Vector<T> v_div(size);
    Vector<T> p_laplace(size);
    size_t idx = 0;
    for (size_t i = 0; i < v_adv.getWidth(); i++) {
        for (size_t j = 0; j < v_adv.getHeight(); i++) {
            for (size_t k = 0; k < v_adv.getDepth(); i++) {
                auto v_div_ijk = v_adv.div(i, j, k);
                auto p_laplace_ijk = p.laplace(i, j, k);
                v_div[idx++] = v_div_ijk;
                p_laplace[idx++] = p_laplace_ijk;
            }
        }
    }
}


int main(int argc, char** argv) {
    // Command line parsing (later)

    // Set up grid
    auto grid = std::make_shared<Grid<ScalarT>>(100, 100, 100, 1);

    // Create velovity and pressure fields
    auto v = std::make_unique<VelocityField<ScalarT>>(grid);
    auto p = std::make_unique<PressureField<ScalarT>>(grid);

    // Init values

    
    double endTime = 10.0;
    double dt = 0.1;

    std::cout << "Initialization done!\n";

    // Time loop
    for (double t = 0; t < endTime; t+=dt) {
        // - Solve advection: implicit euler
        auto v1 = solveAdvection(*v, dt);
        


        // - Pressure corection
        // - Write field (later)
    }

    // - Validation

    // - Visualization

    return 0;
}
