#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, RadiusOfGyration) {
    // Case 1: Valid protein structure
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::pair<bool, double> result = coot::radius_of_gyration_gemmi(st);
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 19.5426, 1e-4);
    }

    // Case 2: Empty structure (should fail)
    {
        gemmi::Structure st;
        std::pair<bool, double> result = coot::radius_of_gyration_gemmi(st);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}