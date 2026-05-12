#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, min_resno_in_chain) {
    // Load standard test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "No model";
    
    // Pick first chain (A)
    ASSERT_FALSE(st.models[0].chains.empty()) << "No chain A";
    const gemmi::Chain& chain = st.models[0].chains[0];
    
    EXPECT_EQ(chain.name, "A");
    
    // Test with valid chain
    std::pair<bool, int> result = coot::minimol::molecule::min_resno_in_chain_gemmi(chain);
    EXPECT_TRUE(result.first);
    EXPECT_EQ(result.second, 1);
    
    // Test with empty chain (edge case)
    gemmi::Chain empty_chain;
    empty_chain.name = "Z";
    std::pair<bool, int> empty_result = coot::minimol::molecule::min_resno_in_chain_gemmi(empty_chain);
    EXPECT_FALSE(empty_result.first);
    EXPECT_EQ(empty_result.second, 99999999);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}