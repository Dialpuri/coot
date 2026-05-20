#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to get chain name for a residue
// Since gemmi::Residue doesn't have a direct chain reference,
// we search through the model's chains to find which one contains the residue
std::string get_chain_name(const gemmi::Structure& st, const gemmi::Residue* res) {
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& r : chain.residues) {
                if (&r == res) {
                    return chain.name;
                }
            }
        }
    }
    return "";
}

TEST(OracleTest, GetResiduePair) {
    // Case 1: Normal case - should find two residues
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        coot::link_by_torsion_t linker;
        std::pair<gemmi::Residue *, gemmi::Residue *> result = 
            linker.get_residue_pair_gemmi(st);
        
        // residue1_found: true
        EXPECT_NE(result.first, nullptr);
        if (result.first) {
            int seqnum1 = result.first->seqid.num.value;
            std::string chain1 = get_chain_name(st, result.first);
            std::string resname1 = result.first->name;
            // residue1: A/1 MET
            EXPECT_EQ(seqnum1, 1);
            EXPECT_EQ(chain1, "A");
            EXPECT_EQ(resname1, "MET");
        }
        
        // residue2_found: true
        EXPECT_NE(result.second, nullptr);
        if (result.second) {
            int seqnum2 = result.second->seqid.num.value;
            std::string chain2 = get_chain_name(st, result.second);
            std::string resname2 = result.second->name;
            // residue2: A/2 GLU
            EXPECT_EQ(seqnum2, 2);
            EXPECT_EQ(chain2, "A");
            EXPECT_EQ(resname2, "GLU");
        }
    }
    
    // Case 2: Structure with no chains - should return null pair
    {
        gemmi::Structure empty_st;
        gemmi::Model model;
        empty_st.models.push_back(model);
        // Note: we don't add any chains to this model
        
        coot::link_by_torsion_t linker;
        std::pair<gemmi::Residue *, gemmi::Residue *> result = 
            linker.get_residue_pair_gemmi(empty_st);
        
        // residue1_found: false
        EXPECT_EQ(result.first, nullptr);
        // residue2_found: false
        EXPECT_EQ(result.second, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}