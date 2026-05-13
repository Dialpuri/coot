#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, sort_residues_by_seqno) {
    // Setup: load the same PDB file as the oracle
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Case 1: chain A residues
    {
        auto& model = st.models[0];
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "A") {
                chain = &c;
                break;
            }
        }
        ASSERT_NE(chain, nullptr) << "chain A not found";
        
        std::vector<gemmi::Residue*> residues;
        for (auto& res : chain->residues) {
            residues.push_back(&res);
        }
        
        std::vector<std::pair<gemmi::Residue*, int>> sorted =
            coot::util::sort_residues_by_seqno(residues.data(), residues.size());
        
        EXPECT_EQ(sorted.size(), 267U) << "chain A should have 267 residues";
        
        int limit = (sorted.size() < 3) ? sorted.size() : 3;
        for (int i = 0; i < limit; ++i) {
            gemmi::Residue* r = sorted[i].first;
            int seqnum = r->seqid.num.value;
            const std::string& resname = r->name;
            if (i == 0) {
                EXPECT_EQ(seqnum, 1);
                EXPECT_EQ(resname, "MET");
            } else if (i == 1) {
                EXPECT_EQ(seqnum, 2);
                EXPECT_EQ(resname, "GLU");
            } else if (i == 2) {
                EXPECT_EQ(seqnum, 3);
                EXPECT_EQ(resname, "ASN");
            }
        }
    }
    
    // Case 2: chain B residues
    {
        auto& model = st.models[0];
        gemmi::Chain* chain = nullptr;
        for (auto& c : model.chains) {
            if (c.name == "B") {
                chain = &c;
                break;
            }
        }
        ASSERT_NE(chain, nullptr) << "chain B not found";
        
        std::vector<gemmi::Residue*> residues;
        for (auto& res : chain->residues) {
            residues.push_back(&res);
        }
        
        std::vector<std::pair<gemmi::Residue*, int>> sorted =
            coot::util::sort_residues_by_seqno(residues.data(), residues.size());
        
        EXPECT_EQ(sorted.size(), 1U) << "chain B should have 1 residue";
        
        int limit = (sorted.size() < 3) ? sorted.size() : 3;
        for (int i = 0; i < limit; ++i) {
            gemmi::Residue* r = sorted[i].first;
            int seqnum = r->seqid.num.value;
            const std::string& resname = r->name;
            if (i == 0) {
                EXPECT_EQ(seqnum, 1);
                EXPECT_EQ(resname, "HOH");
            }
        }
    }
    
    // Case 3: edge case – empty residue list
    {
        std::vector<std::pair<gemmi::Residue*, int>> sorted =
            coot::util::sort_residues_by_seqno(nullptr, 0);
        
        EXPECT_EQ(sorted.size(), 0U) << "empty list should produce empty result";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}