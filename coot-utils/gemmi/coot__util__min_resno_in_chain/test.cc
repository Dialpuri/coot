#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, min_resno_in_chain) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: Valid model - check min resno across all chains
    {
        const gemmi::Model& model = st.models[0];
        auto result = coot::util::min_resno_in_chain_gemmi(&model);
        EXPECT_EQ(result.first, true);
        EXPECT_EQ(result.second, 1);
    }

    // Case 2: NULL model (edge case to test guard clause)
    {
        auto result = coot::util::min_resno_in_chain_gemmi(nullptr);
        EXPECT_EQ(result.first, false);
        EXPECT_EQ(result.second, 99999999);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}