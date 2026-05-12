#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to convert gemmi Residue to identifier for comparison
static coot::residue_identifier get_residue_identifier(const gemmi::Residue& r) {
    return {r.seqid.num.value, r.seqid.icode};
}

TEST(OracleTest, is_a_moving_residue_p) {
    // Load the PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model and chain
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    
    // Need at least 3 residues
    ASSERT_GE(chain.residues.size(), 3u);
    
    gemmi::Residue& res1 = chain.residues[0];
    gemmi::Residue& res2 = chain.residues[1];
    gemmi::Residue& res3 = chain.residues[2];
    
    // Build residues_vec with all residues as moving
    coot::restraints_container_t restraints;
    
    // Add all residues to the moving set
    restraints.residues_vec_moving_set.insert(get_residue_identifier(res1));
    restraints.residues_vec_moving_set.insert(get_residue_identifier(res2));
    restraints.residues_vec_moving_set.insert(get_residue_identifier(res3));
    
    // Case 1: all residues in residues_vec should be in residues_vec_moving_set
    {
        bool r1 = restraints.is_a_moving_residue_p_gemmi(&res1);
        bool r2 = restraints.is_a_moving_residue_p_gemmi(&res2);
        bool r3 = restraints.is_a_moving_residue_p_gemmi(&res3);
        
        EXPECT_TRUE(r1);
        EXPECT_TRUE(r2);
        EXPECT_TRUE(r3);
    }
    
    // Case 2: residue not in residues_vec (should be false)
    if (chain.residues.size() >= 3) {
        gemmi::Residue& res_outside = chain.residues[chain.residues.size() - 1];
        bool r_out = restraints.is_a_moving_residue_p_gemmi(&res_outside);
        EXPECT_FALSE(r_out);
    }
    
    // Case 3: add a residue to residues_vec_moving_set directly and test
    if (chain.residues.size() >= 3) {
        gemmi::Residue& res_outside = chain.residues[chain.residues.size() - 1];
        restraints.residues_vec_moving_set.insert(get_residue_identifier(res_outside));
        bool r_out = restraints.is_a_moving_residue_p_gemmi(&res_outside);
        EXPECT_TRUE(r_out);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}