#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_het_residue_gemmi) {
    // Load the structure using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_FALSE(st.models.empty()) << "No models found";
    
    // Get the first model
    gemmi::Model& model = st.models[0];
    
    // Iterate over chains and residues to find heterogen residues
    bool found_het = false;
    int het_count = 0;
    
    for (gemmi::Chain& chain : model.chains) {
        std::string chain_id = chain.name;
        
        for (gemmi::Residue& res : chain.residues) {
            bool result = coot::molecule_t::is_het_residue_gemmi(&res);
            
            if (result) {
                found_het = true;
                het_count++;
                // Assert the known heterogen residue from oracle output
                EXPECT_EQ(chain_id, "B");
                EXPECT_EQ(res.seqid.num.value, 1);
                EXPECT_EQ(res.name, "HOH");
            }
        }
    }
    
    // Expect exactly 1 heterogen residue (water)
    EXPECT_EQ(het_count, 1);
    EXPECT_TRUE(found_het);
    
    // Edge case: nullptr residue
    bool null_result = coot::molecule_t::is_het_residue_gemmi(nullptr);
    EXPECT_FALSE(null_result);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}