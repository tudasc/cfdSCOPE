#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "Grid.h"
#include "Vector.h"

TEST(TrilerpTests, Basics) {
    Vector<float> init(27, 1.0);

    VelocityField<float> field(std::make_shared<Grid<float>>(3, 3, 3, 1.0), init);

    auto t = field.trilerp({1.5, 1.5, 1.5});
    ASSERT_FLOAT_EQ(t.x, 1.0);
    ASSERT_FLOAT_EQ(t.y, 1.0);
    ASSERT_FLOAT_EQ(t.z, 1.0);   
}