#include <gtest/gtest.h>

#include "Solver.h"
#include "Vector.h"

TEST(SolverTests, ConjugateGradientBasics) {
    std::vector<SparseMatrixEntry<float>> entries = {{0, 0, 1.0}, {1, 1, 1.0}};
    SparseMatrix<float> A(2, 2, entries);
    Vector<float> b = {1.0, 2.0};

    Vector<float> x = pcg(A, b);

    ASSERT_FLOAT_EQ(b[0], x[0]);
    ASSERT_FLOAT_EQ(b[1], x[1]);
}

TEST(SolverTests, ConjugateGradient3x3) {
    std::vector<SparseMatrixEntry<double>> entries = {
        {0, 0, 6.40087},   {0, 1, 0.12948},  {0, 2, 0.171751},
        {1, 0, 0.12948},   {1, 1, 6.78799},  {1, 2, 0.300818},
        {2, 0, 0.1717516}, {2, 1, 0.300818}, {2, 2, 6.17521},
    };
    SparseMatrix<double> A(3, 3, entries);
    Vector<double> b = {0.485769, 0.405699, 0.0174415};

    Vector<double> x = pcg(A, b);

    EXPECT_NEAR(x[0], 0.0747655, 0.001);
    EXPECT_NEAR(x[1], 0.0584342, 0.001);
    EXPECT_NEAR(x[2], -0.00210156, 0.001);
}
 