#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "Matrix.h"
#include "Vector.h"

TEST(LinalgBasicsTest, VectorBasics) {
    Vector<float> a(2);

    a[0] = 1.0;
    a[1] = 3.141;

    ASSERT_EQ(a.getSize(), 2);

    ASSERT_FLOAT_EQ(a[0], 1.0);
    ASSERT_FLOAT_EQ(a[1], 3.141);
}

TEST(LinalgBasicsTest, VectorInitializer) {
    Vector<float> a = {1.0, 2.0, 3.0};

    ASSERT_EQ(a.getSize(), 3);

    ASSERT_FLOAT_EQ(a[0], 1.0);
    ASSERT_FLOAT_EQ(a[1], 2.0);
    ASSERT_FLOAT_EQ(a[2], 3.0);
}

TEST(LinalgBasicsTest, VectorNorm) {
    Vector<float> a = {2.0, 1.0};

    ASSERT_FLOAT_EQ(norm(a), std::sqrt(5.0));
}

TEST(LinalgBasicsTest, VectorDotProduct) {
    Vector<float> a = {1.0, 2.0, 3.0};
    Vector<float> b = {3.0, 2.0, 1.0};

    ASSERT_FLOAT_EQ(dot(a, b), 12.0);
}

TEST(LinalgBasicsTest, MatrixBasics) {
    std::vector<SparseMatrixEntry<float>> entries = {
        {0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 3.0}, {1, 1, 4.0}};

    SparseMatrix<float> m(2, 2, entries);

    ASSERT_EQ(m.getRows(), 2);
    ASSERT_EQ(m.getCols(), 2);

    ASSERT_FLOAT_EQ(m(0, 0), 1.0);
    ASSERT_FLOAT_EQ(m(0, 1), 2.0);
    ASSERT_FLOAT_EQ(m(1, 0), 3.0);
    ASSERT_FLOAT_EQ(m(1, 1), 4.0);
}

TEST(LinalgBasicsTest, MatrixVectorProduct) {
    std::vector<SparseMatrixEntry<float>> entries = {
        {0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 3.0}, {1, 1, 4.0}};

    SparseMatrix<float> m(2, 2, entries);
    Vector<float> a = {1.0, 2.0};
    Vector<float> product = matVecMult(m, a);

    ASSERT_EQ(product.getSize(), m.getRows());
    ASSERT_FLOAT_EQ(product[0], 5.0);
    ASSERT_FLOAT_EQ(product[1], 11.0);
}
