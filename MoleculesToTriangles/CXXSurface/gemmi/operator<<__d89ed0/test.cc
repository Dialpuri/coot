#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, CXXCoord) {
    // Case 1: default-constructed (all zeros)
    {
        CXXCoord<double> c;
        std::ostringstream out;
        stream_CXXCoord_gemmi(out, c);
        EXPECT_EQ(out.str(), "[ 0 0 0 0 ]");
    }

    // Case 2: constructed with x, y, z
    {
        CXXCoord<double> c(1.5, -2.3, 4.7);
        std::ostringstream out;
        stream_CXXCoord_gemmi(out, c);
        EXPECT_EQ(out.str(), "[ 1.5 -2.3 4.7 0 ]");
    }

    // Case 3: constructed with x, y, z, r (all four components)
    {
        CXXCoord<double> c(0.0, 0.0, 0.0, 1.0);
        std::ostringstream out;
        stream_CXXCoord_gemmi(out, c);
        EXPECT_EQ(out.str(), "[ 0 0 0 1 ]");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}