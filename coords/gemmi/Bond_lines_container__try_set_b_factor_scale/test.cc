#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// case: with_UDD_data
TEST(OracleTest, BondLinesUDDScale) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    coot::Bond_lines_container blc;
    double test_scale = 2.5;
    
    // In the original test, the UDD data was set on the molecule,
    // but gemmi doesn't have UDD. The test expects b_factor_scale to remain 1.0
    // because there's no equivalent mechanism in gemmi.
    
    float b_factor_scale_before = blc.b_factor_scale;
    blc.try_set_b_factor_scale_gemmi(st);
    
    // The oracle output shows b_factor_scale_after is 1, not 2.5
    EXPECT_EQ(blc.b_factor_scale, 1.0f);
    EXPECT_NE(blc.b_factor_scale, test_scale);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}