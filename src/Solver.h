#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <limits>

#include "Matrix.h"
#include "Vector.h"

#include "spdlog/spdlog.h"

template <typename T>
Vector<T> jacobi_preconditioner(const SparseMatrix<T>& A, const Vector<T>& x) {
    Vector<T> res(x.getSize());

    for (int i = 0; i < x.getSize(); i++) {
        res[i] = x[i] * (1.0 / A(i, i));
        // res[i] = x[i];
    }

    return res;
}

/**
    Preconditioned conjugate gradient solver
    https://en.wikipedia.org/wiki/Conjugate_gradient_method
*/
template <typename T>
Vector<T> pcg(const SparseMatrix<T>& A, const Vector<T>& b) {
    assert(A.getCols() == A.getRows() &&
           "PCG: Matrix is required to be quadratic.");
    assert(A.getCols() == b.getSize() && "PCG: Dimension mismatch");

    Vector<T> x(b.getSize());
    T tol = 1e-8; // std::numeric_limits<T>::epsilon();

    // Initialize residual vector
    Vector<T> residual = b - A.spmv(x);
    Vector<T> h = jacobi_preconditioner(A, residual);

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
        h = jacobi_preconditioner(A, residual);

        // Update search direction vector
        T newSqrResidNorm = dot(residual, h);
        T beta = newSqrResidNorm / oldSqrResidNorm;
        direction = h + (direction * beta);

        oldSqrResidNorm = newSqrResidNorm;
    }
    SPDLOG_TRACE("PCG Iterations needed: {}", iter);
    SPDLOG_TRACE("Final residual={}", oldSqrResidNorm);
    return x;
}

#endif
