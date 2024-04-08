#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "Grid.h"
#include "Vector.h"

TEST(TrilerpTests, Sanity) {
    Vector<float> init(27*3, 1.0);

    VelocityField<float> field(std::make_shared<Grid<float>>(3, 3, 3, 1.0), init);

    auto t = field.trilerp({1.5, 1.5, 1.5});
    ASSERT_FLOAT_EQ(t.x, 1.0);
    ASSERT_FLOAT_EQ(t.y, 1.0);
    ASSERT_FLOAT_EQ(t.z, 1.0);   
}

TEST(TrilerpTests, SimpleInterpolations) {
    Vector<float> init(27*3, 0.0);

    VelocityField<float> field(std::make_shared<Grid<float>>(3, 3, 3, 1.0), init);

    field.setLeftU(1, 1, 1, 1.0);
    field.setTopV(1, 1, 1, 1.0);
    field.setFrontW(1, 1, 1, 1.0);

    auto t = field.trilerp({1.5, 1.5, 1.5});
    ASSERT_FLOAT_EQ(t.x, 0.5);
    ASSERT_FLOAT_EQ(t.y, 0.5);
    ASSERT_FLOAT_EQ(t.z, 0.5);

    t = field.trilerp({1.0, 1.5, 1.5});
    ASSERT_FLOAT_EQ(t.x, 1.0);
    ASSERT_FLOAT_EQ(t.y, 0.25);
    ASSERT_FLOAT_EQ(t.z, 0.25);

    t = field.trilerp({1.25, 1.25, 1.25});
    ASSERT_FLOAT_EQ(t.x, 0.75*0.75*0.75);
    ASSERT_FLOAT_EQ(t.y, 0.75*0.75*0.75);
    ASSERT_FLOAT_EQ(t.z, 0.75*0.75*0.75);

}