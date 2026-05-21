#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, DeepCopyThisResidue) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    // Find ILE A 10 residue
    gemmi::Residue* target_res = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10 && res.name == "ILE") {
                    target_res = &res;
                    break;
                }
            }
        }
        if (target_res) break;
    }
    ASSERT_NE(target_res, nullptr);

    // Case 1: deep_copy_this_residue without alt_conf (default parameter)
    {
        std::pair<bool, std::string> use_alt_conf = std::make_pair(false, "");
        gemmi::Residue* copied = coot::util::deep_copy_this_residue_gemmi(target_res, use_alt_conf);
        EXPECT_NE(copied, nullptr);
        if (copied) {
            EXPECT_EQ(copied->atoms.size(), 8);
            EXPECT_EQ(copied->seqid.num.value, 10);
            EXPECT_EQ(copied->name, "ILE");
            EXPECT_EQ(copied->seqid.icode, ' ');
            
            // Cleanup the copied residue
            copied->atoms.clear();
            delete copied;
        }
    }
    
    // Case 2: deep_copy_this_residue with alt_conf.first = false (copy all atoms)
    {
        std::pair<bool, std::string> use_alt_conf = std::make_pair(false, "");
        gemmi::Residue* copied = coot::util::deep_copy_this_residue_gemmi(target_res, use_alt_conf);
        EXPECT_NE(copied, nullptr);
        if (copied) {
            EXPECT_EQ(copied->atoms.size(), 8);
            
            copied->atoms.clear();
            delete copied;
        }
    }
    
    // Case 3: deep_copy_this_residue with alt_conf.first = true, empty string (should copy blank altloc atoms)
    {
        std::pair<bool, std::string> use_alt_conf = std::make_pair(true, "");
        gemmi::Residue* copied = coot::util::deep_copy_this_residue_gemmi(target_res, use_alt_conf);
        EXPECT_NE(copied, nullptr);
        if (copied) {
            EXPECT_EQ(copied->atoms.size(), 8);
            
            copied->atoms.clear();
            delete copied;
        }
    }
    
    // Case 4: deep_copy_this_residue with null residue
    {
        std::pair<bool, std::string> use_alt_conf = std::make_pair(false, "");
        gemmi::Residue* copied = coot::util::deep_copy_this_residue_gemmi(nullptr, use_alt_conf);
        EXPECT_EQ(copied, nullptr);
    }
    
    // Case 5: deep_copy_this_residue with alt_conf.first = true and specific altloc
    {
        std::pair<bool, std::string> use_alt_conf = std::make_pair(true, "A");
        gemmi::Residue* copied = coot::util::deep_copy_this_residue_gemmi(target_res, use_alt_conf);
        EXPECT_NE(copied, nullptr);
        if (copied) {
            EXPECT_GT(copied->atoms.size(), 0);
            
            copied->atoms.clear();
            delete copied;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}