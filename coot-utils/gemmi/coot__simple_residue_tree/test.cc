#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SimpleResidueTree) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // case 1: get residue A/10 and build its tree
    {
        gemmi::Residue* res = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& r : chain.residues) {
                    if (chain.name == "A" && r.seqid.num.value == 10) {
                        res = &r;
                        break;
                    }
                }
                if (res) break;
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr) << "Residue not found for chain A/10";
        
        // Get residue info for logging
        gemmi::Chain* chain = nullptr;
        for (auto& model : st.models) {
            for (auto& c : model.chains) {
                for (auto& r : c.residues) {
                    if (&r == res) {
                        chain = &c;
                        break;
                    }
                }
                if (chain) break;
            }
            if (chain) break;
        }
        
        float close_dist_max = 5.0f;
        std::vector<gemmi::Residue*> tree = coot::simple_residue_tree_gemmi(res, st, close_dist_max);
        
        // Verify tree size (from oracle: 268)
        EXPECT_EQ(tree.size(), 268);
        
        // Collect all residues for verification
        std::set<std::pair<int, std::string>> residue_set;
        for (unsigned int i = 0; i < tree.size(); i++) {
            gemmi::Residue* r = tree[i];
            // Find the chain for this residue
            for (auto& model : st.models) {
                for (auto& c : model.chains) {
                    for (auto& rr : c.residues) {
                        if (&rr == r) {
                            residue_set.insert(std::make_pair(r->seqid.num.value, c.name));
                            break;
                        }
                    }
                }
            }
        }
        
        // Check that we have the expected residues from oracle output
        EXPECT_EQ(residue_set.size(), 268);
    }

    // case 2: try residue A/20
    {
        gemmi::Residue* res = nullptr;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& r : chain.residues) {
                    if (chain.name == "A" && r.seqid.num.value == 20) {
                        res = &r;
                        break;
                    }
                }
                if (res) break;
            }
            if (res) break;
        }
        ASSERT_NE(res, nullptr) << "Residue not found for chain A/20";
        
        float close_dist_max = 4.0f;
        std::vector<gemmi::Residue*> tree = coot::simple_residue_tree_gemmi(res, st, close_dist_max);
        
        // Verify tree size (from oracle: 268)
        EXPECT_EQ(tree.size(), 268);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}