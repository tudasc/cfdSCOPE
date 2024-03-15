#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

#include "Vector.h"

/**
    Data type representing one entry in the sparse matrix
*/
template <typename T>
class SparseMatrixEntry {
    size_t row;
    size_t col;
    T value;
};

/**
  Sparse matrix in compressed-sparse-row format (CSR)

  c.f., https://en.wikipedia.org/wiki/Sparse_matrix
*/
template <typename T>
class SparseMatrix {

  public:
    /**
        Constructor to build sparse matrix from a vector of entries
    */
    SparseMatrix(size_t rows, size_t cols,
                 std::vector<SparseMatrixEntry<T>> entries)
        : _rows(rows), _cols(cols), _v(entries.size()),
          _colIndex(entries.size()), _rowIndex(rows + 1) {

        // sort matrix entries for easier construction
        std::sort(entries.begin(), entries.end(),
                  [](SparseMatrixEntry<T> a, SparseMatrixEntry<T> b) {
                      if (a.row != b.row)
                          return a.row > b.row;
                      else
                          return a.col > b.col;
                  });

        // iterate over entries and put it in the data structure by adjusting the CSR-arrays
        size_t lastRow = 0;
        for (size_t i = 0; i < entries.size(); i++) {
            SparseMatrixEntry<T>& entry = entries[i];

            assert(entry.value != 0.0 && "By definition, a sparse matrix only "
                                         "contains non-zero elements.");
            assert((entry.row >= 0 && entry.row < _rows && entry.col >= _rows &&
                    entry.col < _cols) &&
                   "Out-of-bounds sparse-matrix element!");

            _v[i] = entry.value;
            _colIndex[i] = entry.col;

            // are we in a new row?
            if (entry.row > lastRow) {
                // set row index of new row (and all possibly skipped rows) to
                // the current index
                for (size_t j = lastRow; j <= entry.row; j++) {
                    _rowIndex[j] = i;
                }

                lastRow = entry.row;
            }
        }

        _rowIndex[rows] = entries.size();
    }

    /**
        Matrix access operator
    */
    const T& operator()(size_t r, size_t c) {
        assert((r >= 0 && r < _rows && c >= 0 && c < _cols) &&
               "Out-of-bounds access!");

        size_t row_start = _rowIndex[r];
        size_t row_end = _rowIndex[r + 1];

        for (int i = row_start; i < row_end; i++) {
            if (_colIndex[i] < c)
                continue;
            if (_colIndex[i] == c)
                return _v[i];
            if (_colIndex[i] > c)
                return 0.0;
        }
        return 0.0;
    }

    size_t getCols() const { return _cols; }

    size_t getRows() const { return _rows; }

  private:
    // dimensions
    size_t _rows;
    size_t _cols;

    // CSR arrays
    Vector<T> _v;
    Vector<size_t> _colIndex;
    Vector<size_t> _rowIndex;
};

/**
    Matrix-vector product
*/
template <typename T>
SparseMatrix<T> matVecMult(const SparseMatrix<T>& A, const Vector<T>& v) {
    assert(A.getCols() == v.getSize() && "matVecMult: Dimension mismatch!");

    Vector<T> res(A.getRows());

    for (size_t i = 0; i < A.getRows(); i++) {
        for (size_t j = 0; j < A.getCols(); j++) {
            res[i] += A(i, j) * v[j];
        }
    }

    return res;
}
