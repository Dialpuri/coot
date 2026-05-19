#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to get a residue from the structure using chain ID and sequence number
static gemmi::Residue* get_residue_by_seqid(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, DeleteTheMatchedResiduesMatchedResidue) {
    // Test case 1: delete residue_1 (short_fragment_is_in_first_selection = false)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Get two residues from chain A
        gemmi::Residue* res1 = get_residue_by_seqid(st.models[0], "A", 10);
        gemmi::Residue* res2 = get_residue_by_seqid(st.models[0], "A", 20);
        
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        int nres_before = 0;
        for (const auto& chain : st.models[0].chains) {
            nres_before += chain.residues.size();
        }
        EXPECT_EQ(nres_before, 268);
        
        coot::delete_the_matched_residues_matched_residue_gemmi(st, res1, res2, false);
        
        int nres_after = 0;
        for (const auto& chain : st.models[0].chains) {
            nres_after += chain.residues.size();
        }
        EXPECT_EQ(nres_after, 267);
    }
    
    // Test case 2: delete residue_2 (short_fragment_is_in_first_selection = true)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        // Get two residues from chain A
        gemmi::Residue* res1 = get_residue_by_seqid(st.models[0], "A", 10);
        gemmi::Residue* res2 = get_residue_by_seqid(st.models[0], "A", 20);
        
        ASSERT_NE(res1, nullptr);
        ASSERT_NE(res2, nullptr);
        
        int nres_before = 0;
        for (const auto& chain : st.models[0].chains) {
            nres_before += chain.residues.size();
        }
        EXPECT_EQ(nres_before, 268);
        
        coot::delete_the_matched_residues_matched_residue_gemmi(st, res1, res2, true);
        
        int nres_after = 0;
        for (const auto& chain : st.models[0].chains) {
            nres_after += chain.residues.size();
        }
        EXPECT_EQ(nres_after, 267);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}