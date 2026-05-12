#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeepCopyResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: valid residue
    {
        // Find residue with chain A, sequence number 10
        gemmi::Residue *original = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (chain.name == "A" && res.seqid.num.value == 10) {
                        original = &res;
                        break;
                    }
                }
                if (original) break;
            }
            if (original) break;
        }
        
        EXPECT_TRUE(original != nullptr) << "original residue should be found";
        
        if (original) {
            size_t original_atoms = original->atoms.size();
            gemmi::Residue *copy = coot::util::deep_copy_this_residue_gemmi(original);
            EXPECT_TRUE(copy != nullptr) << "copy should be created";
            
            if (copy) {
                size_t copy_atoms = copy->atoms.size();
                EXPECT_EQ(original->seqid.num.value, 10) << "original_seqnum";
                EXPECT_EQ(copy->seqid.num.value, 10) << "copy_seqnum";
                EXPECT_EQ(original->seqid.icode, copy->seqid.icode) << "inscode";
                EXPECT_EQ(original->name, copy->name) << "resname";
                EXPECT_EQ(original_atoms, copy_atoms) << "atoms_match";
                EXPECT_EQ(original->seqid.num.value, copy->seqid.num.value) << "seqnum_match";
                EXPECT_EQ(original->seqid.icode, copy->seqid.icode) << "inscode_match";
                EXPECT_EQ(original->name, copy->name) << "resname_match";
                delete copy;
            }
        }
    }

    // case 2: NULL input — should return NULL
    {
        gemmi::Residue *null_res = nullptr;
        gemmi::Residue *copy = coot::util::deep_copy_this_residue_gemmi(null_res);
        EXPECT_TRUE(copy == nullptr) << "copy_from_null should be nullptr";
    }

    // case 3: non-existent residue — should return NULL
    {
        gemmi::Residue *original = nullptr;
        // Try to find residue with sequence number 9999 (should not exist)
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (chain.name == "A" && res.seqid.num.value == 9999) {
                        original = &res;
                        break;
                    }
                }
                if (original) break;
            }
            if (original) break;
        }
        
        gemmi::Residue *copy = coot::util::deep_copy_this_residue_gemmi(original);
        EXPECT_TRUE(original == nullptr) << "original_from_invalid_cid should be nullptr";
        EXPECT_TRUE(copy == nullptr) << "copy_from_invalid should be nullptr";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}