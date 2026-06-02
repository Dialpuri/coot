#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, Cartesian) {
    // Case 1: typical Cartesian with all four components
    {
        Cartesian_gemmi c(1.5, -2.3, 4.7, 0.85);
        EXPECT_NEAR(c.get_x(), 1.5, 1e-4);
        EXPECT_NEAR(c.get_y(), -2.3, 1e-4);
        EXPECT_NEAR(c.get_z(), 4.7, 1e-4);
        EXPECT_NEAR(c.get_a(), 0.85, 1e-4);

        std::ostringstream oss;
        operator_stream_gemmi(oss, c);
        EXPECT_EQ(oss.str(), "  1.500000  -2.300000   4.700000   0.850000");
    }

    // Case 2: Cartesian with default 'a' (1.0)
    {
        Cartesian_gemmi c(0.0, 0.0, 0.0);
        EXPECT_NEAR(c.get_x(), 0.0, 1e-4);
        EXPECT_NEAR(c.get_y(), 0.0, 1e-4);
        EXPECT_NEAR(c.get_z(), 0.0, 1e-4);
        EXPECT_NEAR(c.get_a(), 1.0, 1e-4);

        std::ostringstream oss;
        operator_stream_gemmi(oss, c);
        EXPECT_EQ(oss.str(), "  0.000000   0.000000   0.000000   1.000000");
    }

    // Case 3: edge case — large and negative values
    {
        Cartesian_gemmi c(-999.123456, 12345.0, -0.000001, 2.0);
        EXPECT_NEAR(c.get_x(), -999.123456, 1e-4);
        EXPECT_NEAR(c.get_y(), 12345.0, 1e-4);
        EXPECT_NEAR(c.get_z(), -0.000001, 1e-4);
        EXPECT_NEAR(c.get_a(), 2.0, 1e-4);

        std::ostringstream oss;
        operator_stream_gemmi(oss, c);
        EXPECT_EQ(oss.str(), "-999.123456 12345.000000  -0.000001   2.000000");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}