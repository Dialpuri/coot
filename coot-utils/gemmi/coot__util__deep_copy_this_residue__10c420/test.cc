#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeepCopyResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue - copy MET 1 in chain A
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        
        EXPECT_NE(chain_a, nullptr) << "Chain A should exist";
        
        const gemmi::Residue* original = nullptr;
        if (chain_a) {
            for (const gemmi::Residue& res : chain_a->residues) {
                if (res.seqid.num.value == 1) {
                    original = &res;
                    break;
                }
            }
        }
        
        EXPECT_NE(original, nullptr) << "Original residue should not be null";
        
        if (original) {
            // Before copy - get original info
            int seqNum_orig = original->seqid.num.value;
            std::string resName_orig = original->name;
            
            int nAtoms_orig = static_cast<int>(original->atoms.size());
            
            // Make the copy
            gemmi::Residue* copy = coot::util::deep_copy_this_residue_gemmi(original);
            
            EXPECT_NE(copy, nullptr) << "Copy should be created";
            
            if (copy) {
                // After copy - get copy info
                int seqNum_copy = copy->seqid.num.value;
                std::string resName_copy = copy->name;
                
                int nAtoms_copy = static_cast<int>(copy->atoms.size());
                
                EXPECT_EQ(seqNum_orig, seqNum_copy);
                EXPECT_EQ(resName_orig, resName_copy);
                EXPECT_EQ(nAtoms_orig, nAtoms_copy);
                
                // Verify properties match
                bool same_seqnum = (seqNum_orig == seqNum_copy);
                bool same_resname = (resName_orig == resName_copy);
                bool same_natoms = (nAtoms_orig == nAtoms_copy);
                
                EXPECT_TRUE(same_seqnum && same_resname && same_natoms);
                
                // Clean up
                delete copy;
            }
        }
    }

    // case 2: NULL input - should return NULL
    {
        const gemmi::Residue* null_res = nullptr;
        gemmi::Residue* copy = coot::util::deep_copy_this_residue_gemmi(null_res);
        EXPECT_EQ(copy, nullptr) << "Copy should be null for null input";
    }

    // case 3: another valid residue - GLU 5 in chain A
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain* chain_a = nullptr;
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
        }
        
        EXPECT_NE(chain_a, nullptr) << "Chain A should exist";
        
        const gemmi::Residue* original = nullptr;
        if (chain_a) {
            for (const gemmi::Residue& res : chain_a->residues) {
                if (res.seqid.num.value == 5) {
                    original = &res;
                    break;
                }
            }
        }
        
        EXPECT_NE(original, nullptr) << "Original residue should not be null";
        
        if (original) {
            int nAtoms_orig = static_cast<int>(original->atoms.size());
            
            gemmi::Residue* copy = coot::util::deep_copy_this_residue_gemmi(original);
            
            EXPECT_NE(copy, nullptr) << "Copy should be created";
            
            if (copy) {
                int nAtoms_copy = static_cast<int>(copy->atoms.size());
                
                EXPECT_EQ(5, copy->seqid.num.value);
                EXPECT_STREQ("GLN", copy->name.c_str());
                EXPECT_EQ(nAtoms_orig, nAtoms_copy);
                EXPECT_TRUE(nAtoms_orig == nAtoms_copy);
                
                delete copy;
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}