#include <cstddef>

#include "Vector.h"

template<typename T>
class SparseMatrix {
size_t rows;
size_t cols;




public:
    SparseMatrix(size_t rows, size_t cols);

    T& operator()(size_t r, size_t c) {

    }

    size_t getCols() const {
        return cols;
    }

    size_t getRows() const {
        return rows;
    }



};

template<typename T>
SparseMatrix<T> matVecMult(const SparseMatrix<T>& A, const Vector<T>& v ) {

}

template<typename T>
T dot(const Vector<T>& u, const Vector<T>& v ) {

}

