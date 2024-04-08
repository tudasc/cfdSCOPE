#ifndef UTIL_H
#define UTIL_H

#include "Grid.h"
#include "Matrix.h"
#include "Vector.h"

#include <iostream>

template <typename T>
inline void dumpMatrix(const SparseMatrix<T>& A) {
    // std::cout << "[";
    for (auto i = 0; i < A.getRows(); i++) {
        // std::cout << "[";
        for (auto j = 0; j < A.getCols(); j++) {
            std::cout << A(i, j);
            if (j < A.getCols() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "\n";
        if (i < A.getRows() - 1) {
            // std::cout << "]\n";
        } else {

            // std::cout << "]]\n";
        }
    }
}

template <typename T>
inline void dumpVector(const Vector<T>& v) {
    std::cout << "[ ";
    for (auto i = 0; i < v.getSize(); i++) {
        std::cout << v[i];
        if (i < v.getSize() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << " ]\n";
}

template <typename T>
inline void dumpVectorComponent(const Vector<T>& v, int comp, int stride) {
    std::cout << "[ ";
    for (auto i = comp; i < v.getSize(); i += stride) {
        std::cout << v[i];
        if (i < v.getSize() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << " ]\n";
}

template <typename T>
inline void dumpField(const VelocityField<T>& v) {}

#endif // UTIL_H
