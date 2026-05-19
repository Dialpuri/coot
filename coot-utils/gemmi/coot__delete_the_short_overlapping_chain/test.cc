#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeleteTheShortOverlappingChain) {
    using namespace coot;

    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test case 1: DELETE_FIRST_CHAIN (chain A)
    {
        std::cout << "=== Test 1: DELETE_FIRST_CHAIN (chain A) ===" << std::endl;
        
        // Before: count chains
        int n_chains_before = st.models[0].chains.size();
        EXPECT_EQ(n_chains_before, 2);
        
        // Build delete_a_chain_t with DELETE_FIRST_CHAIN
        delete_a_chain_t dac(true, true, true);
        dac.delete_type = delete_a_chain_t::DELETE_FIRST_CHAIN;
        
        // Call the function
        delete_the_short_overlapping_chain_gemmi(dac, st, "A", "B");
        
        // After: count chains
        int n_chains_after = st.models[0].chains.size();
        EXPECT_EQ(n_chains_after, 1);
        
        // Check if chain A was deleted
        bool chain_a_exists = false;
        bool chain_b_exists = false;
        for (const auto& chain : st.models[0].chains) {
            if (chain.name == "A") chain_a_exists = true;
            if (chain.name == "B") chain_b_exists = true;
        }
        EXPECT_FALSE(chain_a_exists);
        EXPECT_TRUE(chain_b_exists);
    }
    
    // Reload structure for second test
    st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test case 2: DELETE_SECOND_CHAIN (chain B)
    {
        std::cout << "=== Test 2: DELETE_SECOND_CHAIN (chain B) ===" << std::endl;
        
        // Before: count chains
        int n_chains_before = st.models[0].chains.size();
        EXPECT_EQ(n_chains_before, 2);
        
        // Build delete_a_chain_t with DELETE_SECOND_CHAIN
        delete_a_chain_t dac(true, true, true);
        dac.delete_type = delete_a_chain_t::DELETE_SECOND_CHAIN;
        
        // Call the function
        delete_the_short_overlapping_chain_gemmi(dac, st, "A", "B");
        
        // After: count chains
        int n_chains_after = st.models[0].chains.size();
        EXPECT_EQ(n_chains_after, 1);
        
        // Check if chain B was deleted
        bool chain_a_exists = false;
        bool chain_b_exists = false;
        for (const auto& chain : st.models[0].chains) {
            if (chain.name == "A") chain_a_exists = true;
            if (chain.name == "B") chain_b_exists = true;
        }
        EXPECT_TRUE(chain_a_exists);
        EXPECT_FALSE(chain_b_exists);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}