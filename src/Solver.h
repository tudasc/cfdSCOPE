#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <limits>

#include "Matrix.h"
#include "Vector.h"

#include "spdlog/spdlog.h"

/**
    Conjugate gradient solver
    https://en.wikipedia.org/wiki/Conjugate_gradient_method
*/
template <typename T>
Vector<T> pcg(const SparseMatrix<T>& A, const Vector<T>& b) {
    Vector<T> x(b.getSize());
    T tol = 1e-8; // std::numeric_limits<T>::epsilon();

    auto maxIter = 100;

    // Initialize residual vector
    Vector<T> residual = b - A.spmv(x);
    // Initialize search direction vector
    Vector<T> direction = residual;

    T oldSqrResidNorm = dot(residual, residual);

    // Iterate until convergence
    while (oldSqrResidNorm > tol) {
        if (maxIter-- <= 0) {
            spdlog::error("Max iterations reached - aborting...");
            assert(false && "PCG did not converge");
        }
        SPDLOG_TRACE("PCG residual={}", oldSqrResidNorm);
        Vector<T> z = A.spmv(direction);

        T step_size = dot(residual, residual) / dot(direction, z);

        // Update solution
        x = x + (direction * step_size);
        // Update residual
        residual = residual - (z * step_size);
        T newSqrResidNorm = dot(residual, residual);

        // Update search direction vector
        T beta = newSqrResidNorm / oldSqrResidNorm;
        direction = residual + (direction * beta);

        oldSqrResidNorm = newSqrResidNorm;
    }
    SPDLOG_TRACE("Final residual={}", oldSqrResidNorm);
    return x;
}

#endif
