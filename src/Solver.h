#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <limits>

#include "Matrix.h"
#include "Vector.h"

#include "spdlog/spdlog.h"

template<typename T>
struct Preconditioner {
    virtual Vector<T> precondition(const SparseMatrix<T>& A, const Vector<T>& x) {
        return x;
    }
};

template <typename T>
struct JacobiPreconditioner: public Preconditioner<T> {
    
    JacobiPreconditioner(const SparseMatrix<T>& A) : _rDiag(A.getRows()) {
        assert(A.getCols() == A.getRows());

        auto n = _rDiag.getSize();
        for (int i = 0; i < n; i++) {
            _rDiag[i] = 1.0 / A(i, i);
        }
    }

    Vector<T> precondition(const SparseMatrix<T>& A, const Vector<T>& x) {
        auto n = _rDiag.getSize();
        assert(x.getSize() == n);

        Vector<T> res(n);
        for (int i = 0; i < n; i++) {
            res[i] = _rDiag[i] * x[i];
        }
        return res;
    }

private:
    Vector<T> _rDiag;

};

template <typename T>
struct DICPreconditioner: public Preconditioner<T> {
    
    /**
    * Based on OpenFOAM's simplified diagonal-based incomplete Cholesky preconditioner: 
    * https://develop.openfoam.com/Development/openfoam/-/blob/OpenFOAM-v2112/src/OpenFOAM/matrices/lduMatrix/preconditioners/DICPreconditioner/DICPreconditioner.C
    */
    DICPreconditioner(const SparseMatrix<T>& A) : _rDiag(A.getRows()) {
        assert(A.getCols() == A.getRows());

        auto n = _rDiag.getSize();
        for (int i = 0; i < n; i++) {
            _rDiag[i] = A(i, i);
        }

        // Calculate the DIC diagonal
        for (size_t i = 0; i < n; i++) {
            for (size_t idx = A.getRowIndex(i); idx < A.getRowIndex(i + 1); idx++) {
                size_t j = A.getColIndex(idx);
                if (j > i) {
                    _rDiag[j] -= A.getVal(idx) * A.getVal(idx) / _rDiag[i];
                }
            }
            // for (int j = i+1; j < n; j++) {
            //     auto a_ij = A(i,j);
            //     if (a_ij != 0) {
            //         _rDiag[j] -= a_ij * a_ij / _rDiag[i];
            //     }
            // }
        }

        // Calculate the reciprocal of the preconditioned diagonal
        for (int i = 0; i < n; i++) {
            _rDiag[i] = 1.0 / _rDiag[i];
        }
    }

    Vector<T> precondition(const SparseMatrix<T>& A, const Vector<T>& x) {
        auto n = _rDiag.getSize();
        assert(x.getSize() == n);

        Vector<T> res(n);
        for (int i = 0; i < n; i++) {
            res[i] = _rDiag[i] * x[i];
        }

        // Apply upper part of preconditioner
        for (int i = 0; i < n; i++) {
            for (int idx = A.getRowIndex(i); idx < A.getRowIndex(i + 1); idx++) {
                size_t j = A.getColIndex(idx);
                if (j > i) {
                    res[j] -= _rDiag[j] * A.getVal(idx) * res[i];
                }
            }
            // for (int j = i+1; j < n; j++) {
            //     auto a_ij = A(i,j);
            //     if (a_ij != 0) {
            //         res[j] -= _rDiag[j] * a_ij * res[i];
            //     }
            // }   
        }

        // Apply lower part of preconditioner
        for (int i = n-1; i >= 0; i--) {
            // spdlog::warn("i: {}, a: {}, b: {}", i, A.getRowIndex(i + 1)-1, A.getRowIndex(i));
            for (int idx = A.getRowIndex(i + 1) - 1; idx >= (int) A.getRowIndex(i); idx--) {
                // spdlog::warn("idx: {}, ri: {}", idx, A.getRowIndex(i));
                size_t j = A.getColIndex(idx);
                if (j > i) {
                    res[i] -= _rDiag[i] * A.getVal(idx) * res[j];
                }
            }
            // for (int j = n-1; j >= i+1; j--) {
            //     auto a_ij = A(i,j);
            //     if (a_ij != 0) {
            //         res[i] -= _rDiag[i] * a_ij * res[j];
            //     }
            // }   
        }
        return res;
    }

private:
    Vector<T> _rDiag;

};


// template <typename T>
// Vector<T> dic_preconditioner(const SparseMatrix<T>& A, const Vector<T>& x) {
//     Vector<T> res(x.getSize());

//     Vector<T> rDiag(x.getSize());

//     for (int i = 0; i < x.getSize(); i++) {
//         rDiag[i] = A(i, i);
//     }

//     for (int i = 0; i < x.getSize(); i++) {
//         for (int j = i+1; j < x.getSize(); j++) {
//             auto a_ij = A(i,j);
//             if (a_ij != 0) {
//                 rDiag[j] -= a_ij * a_ij / rDiag[i];
//             }
//         }
//     }

//     for (int i = 0; i < x.getSize(); i++) {
//         rDiag[i] = 1.0 / rDiag[i];
//     }

//     // Apply preconditioner
//     for (int i = 0; i < x.getSize(); i++) {
//         res[i] = rDiag[i] * x[i];
//     }

//     for (int i = 0; i < x.getSize(); i++) {
//         for (int j = i+1; j < x.getSize(); j++) {
//             auto a_ij = A(i,j);
//             if (a_ij != 0) {
//                 res[j] -= rDiag[j] * a_ij * res[i];
//             }
//         }   
//     }

//      for (int i = x.getSize()-1; i >= 0; i--) {
//         for (int j = x.getSize() -1; j >= i+1; j--) {
//             auto a_ij = A(i,j);
//             if (a_ij != 0) {
//                 res[i] -= rDiag[i] * a_ij * res[j];
//             }
//         }   
//     }


//     return res;
// }

/**
    Preconditioned conjugate gradient solver
    https://en.wikipedia.org/wiki/Conjugate_gradient_method
*/
template <typename T>
Vector<T> pcg(const SparseMatrix<T>& A, const Vector<T>& b, Preconditioner<T>& p, const Vector<T>& guess = {}) {
    assert(A.getCols() == A.getRows() &&
           "PCG: Matrix is required to be quadratic.");
    assert(A.getCols() == b.getSize() && "PCG: Dimension mismatch");

    if (!A.isSymmetric()) {
        spdlog::warn("PCG: Matrix is not symmetric!");
    }

    Vector<T> x = Vector<T>(b.getSize());
    if (guess.getSize() == b.getSize()) {
        x = guess;
    }

    T tol = 1e-8; // std::numeric_limits<T>::epsilon();

    // Initialize residual vector
    Vector<T> residual = b - A.spmv(x);
    Vector<T> h = p.precondition(A, residual);

    // Initialize search direction vector
    Vector<T> direction = h;

    T oldSqrResidNorm = dot(residual, h);

    unsigned int iter = 0;

    // Iterate until convergence
    while (dot(residual, residual) > tol) {
        if (iter++ > 1000) {
            spdlog::error("Max iterations reached - aborting...");
            assert(false && "PCG did not converge");
        }
        SPDLOG_TRACE("PCG residual={}", oldSqrResidNorm);
        Vector<T> z = A.spmv(direction);

        T step_size = dot(residual, h) / dot(direction, z);

        // Update solution
        x = x + (direction * step_size);
        // Update residual
        residual = residual - (z * step_size);
        h = p.precondition(A, residual);

        // Update search direction vector
        T newSqrResidNorm = dot(residual, h);
        T beta = newSqrResidNorm / oldSqrResidNorm;
        direction = h + (direction * beta);

        oldSqrResidNorm = newSqrResidNorm;
    }
    spdlog::debug("PCG Iterations needed: {}", iter);
    SPDLOG_TRACE("Final residual={}", oldSqrResidNorm);
    return x;
}

#endif
