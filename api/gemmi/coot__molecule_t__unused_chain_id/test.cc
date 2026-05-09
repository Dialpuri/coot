#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, UnusedChainId) {
    // Load the test molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    const gemmi::Model& model = st.models[0];
    
    // Call unused_chain_id_gemmi()
    std::pair<bool, std::string> result = coot::molecule_t::unused_chain_id_gemmi(model);
    
    EXPECT_EQ(result.first, true);
    EXPECT_EQ(result.second, "C");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}