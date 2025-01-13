#include "IO.h"
#include <cstdio>
#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(IoBasicTest, BasicAssertions) {

    std::string filename = std::string(TEST_RESOURCE_DIR) + "/sample_io";
    auto pair = read_from_file<double>(filename);
    auto filename_out = filename + "_out";
    write_to_raw_file(pair.first, pair.second, filename_out);

    std::ifstream f1(filename_out);
    std::ifstream f2(filename);

    ASSERT_TRUE(std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                           std::istreambuf_iterator<char>(),
                           std::istreambuf_iterator<char>(f2.rdbuf())));
    // clean up
    std::remove(filename_out.c_str());
}
