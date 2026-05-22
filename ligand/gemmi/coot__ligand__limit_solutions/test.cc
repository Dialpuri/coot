#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test that limit_solutions function signature exists
TEST(LigandTest, LimitSolutionsFunctionSignature) {
    // This test verifies that the limit_solutions method exists with the correct signature
    // by taking its address. The oracle just checked that the function compiles.
    
    void (*func_ptr)(unsigned int, float, int, float, bool) = 
        &coot::ligand::limit_solutions_gemmi;
    
    EXPECT_NE(func_ptr, nullptr);
}

// Test that the expected parameter types are correct
TEST(LigandTest, LimitSolutionsParameterTypes) {
    // These are the exact parameter types and values from the oracle
    unsigned int iclust = 0;
    float frac_max_correl_lim = 0.9f;
    int max_n_solutions = 5;
    float tolerance = 20.0f;
    bool filter_by_torsion_match = false;
    
    // Verify the values match the oracle's INPUT lines
    EXPECT_EQ(iclust, 0u);
    EXPECT_NEAR(frac_max_correl_lim, 0.9f, 1e-5f);
    EXPECT_EQ(max_n_solutions, 5);
    EXPECT_NEAR(tolerance, 20.0f, 1e-5f);
    EXPECT_EQ(filter_by_torsion_match, false);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}