#include <gtest/gtest.h>
#include <sstream>
#include "function.hh"

TEST(OracleTest, CXXCoordOutput) {
    // case: default constructed CXXCoord (all zeros)
    {
        CXXCoord<double> coord1;
        std::ostringstream oss1;
        oss1 << coord1;
        EXPECT_EQ(oss1.str(), "[ 0 0 0 0 ]");
    }

    // case: 3-argument constructor (x, y, z)
    {
        CXXCoord<double> coord2(1.0, 2.0, 3.0);
        std::ostringstream oss2;
        oss2 << coord2;
        EXPECT_EQ(oss2.str(), "[ 1 2 3 0 ]");
    }

    // case: 4-argument constructor (x, y, z, r)
    {
        CXXCoord<float> coord3(1.5f, 2.5f, 3.5f, 4.5f);
        std::ostringstream oss3;
        oss3 << coord3;
        EXPECT_EQ(oss3.str(), "[ 1.5 2.5 3.5 4.5 ]");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}