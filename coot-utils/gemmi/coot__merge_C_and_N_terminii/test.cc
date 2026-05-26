#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "function.hh"

// Test the merge_C_and_N_terminii function with various parameters
TEST(OracleTest, mergeCAndNTermini) {
    // Setup: load PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";
    
    // Create dummy xmap
    clipper::Spacegroup spgr;
    clipper::Cell cell;
    clipper::Xmap<float> xmap;
    
    // Test 1: call with default parameters (no symmetry, no missing loop fit)
    {
        EXPECT_NO_THROW(coot::merge_C_and_N_terminii_gemmi(st, xmap, false, false)) 
            << "merge_C_and_N_terminii should complete without exceptions with use_symmetry=false, using_missing_loop_fit=false";
    }
    
    // Test 2: call with symmetry enabled
    {
        EXPECT_NO_THROW(coot::merge_C_and_N_terminii_gemmi(st, xmap, true, false)) 
            << "merge_C_and_N_terminii should complete without exceptions with use_symmetry=true, using_missing_loop_fit=false";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}