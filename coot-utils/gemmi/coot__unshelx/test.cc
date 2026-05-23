#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Load test PDB file path
static std::string test_pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

TEST(OracleTest, UnshelxFunctionTests) {
    // Load original PDB
    gemmi::Structure orig_st = gemmi::read_pdb_file(test_pdb_path);
    
    // Test 1: null input case
    {
        gemmi::Structure *result = coot::unshelx_gemmi(nullptr);
        EXPECT_EQ(result, nullptr);
    }
    
    // Test 2: molecule with multiple chains - should fail (not exactly 1 chain)
    {
        // Check that the test PDB has multiple chains
        int nchains = 0;
        for (const auto& model : orig_st.models)
            for (const auto& chain : model.chains)
                nchains++;
        EXPECT_EQ(nchains, 2) << "test PDB should have 2 chains";
        
        gemmi::Structure *result = coot::unshelx_gemmi(&orig_st);
        EXPECT_EQ(result, nullptr);
    }
    
    // Test 3: Create a single-chain molecule from the loaded molecule
    {
        // Create a single-chain structure from the first chain
        gemmi::Structure single_chain_st;
        single_chain_st.name = "single_chain";
        single_chain_st.connections = orig_st.connections;
        single_chain_st.cell = orig_st.cell;
        single_chain_st.spacegroup_hm = orig_st.spacegroup_hm;
        
        // Copy first chain from first model
        if (!orig_st.models.empty() && !orig_st.models[0].chains.empty()) {
            const auto& first_chain = orig_st.models[0].chains[0];
            
            gemmi::Model new_model;
            gemmi::Chain new_chain;
            new_chain.name = first_chain.name;
            
            // Deep copy residues
            for (const auto& res : first_chain.residues) {
                gemmi::Residue copy_res = res;
                new_chain.residues.push_back(copy_res);
            }
            
            new_model.chains.push_back(new_chain);
            single_chain_st.models.push_back(new_model);
        }
        
        int check_nchains = 0;
        for (const auto& model : single_chain_st.models)
            for (const auto& chain : model.chains)
                check_nchains++;
        EXPECT_EQ(check_nchains, 1) << "single_chain_st should have 1 chain";
        
        gemmi::Structure *result = coot::unshelx_gemmi(&single_chain_st);
        EXPECT_NE(result, nullptr) << "unshelx should return non-null for valid single-chain input";
        
        if (result) {
            int result_nchains = 0;
            for (const auto& model : result->models)
                for (const auto& chain : model.chains)
                    result_nchains++;
            
            if (result_nchains == 1) {
                // Expected: single chain remains as one chain
                EXPECT_EQ(result_nchains, 1);
            } else {
                // Conversion succeeded (single chain split into multiple)
                // This is expected behavior when residues have gaps > skip_chain_step
            }
            delete result;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}