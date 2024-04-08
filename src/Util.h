#ifndef UTIL_H
#define UTIL_H

#include "Grid.h"
#include "Matrix.h"
#include "Vector.h"

#include <iostream>
#include <sstream>

template <typename T>
inline std::string dumpMatrix(const SparseMatrix<T>& A) {
    std::stringstream s;
    // std::cout << "[";
    for (auto i = 0; i < A.getRows(); i++) {
        // std::cout << "[";
        for (auto j = 0; j < A.getCols(); j++) {
            s << A(i, j);
            if (j < A.getCols() - 1) {
                s << ", ";
            }
        }
        s << "\n";
        if (i < A.getRows() - 1) {
            // std::cout << "]\n";
        } else {

            // std::cout << "]]\n";
        }
    }
    return s.str();
}

template <typename T>
inline std::string dumpVector(const Vector<T>& v) {
    std::stringstream s;
    s << "[ ";
    for (auto i = 0; i < v.getSize(); i++) {
        s << v[i];
        if (i < v.getSize() - 1) {
            s << ", ";
        }
    }
    s << " ]\n";
    return s.str();
}

template <typename T>
inline std::string dumpVectorComponent(const Vector<T>& v, int comp,
                                       int stride) {
    std::stringstream s;
    s << "[ ";
    for (auto i = comp; i < v.getSize(); i += stride) {
        s << v[i];
        if (i < v.getSize() - 1) {
            s << ", ";
        }
    }
    s << " ]\n";
    return s.str();
}

#endif // UTIL_H
