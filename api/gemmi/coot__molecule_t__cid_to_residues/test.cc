#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cid_to_residues_cases) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: Valid single-residue CID (should find 1 residue)
    {
        std::string cid = "//A/10";
        std::vector<gemmi::Residue *> residues = coot::cid_to_residues_gemmi(cid, st);
        EXPECT_EQ(residues.size(), 1u);
        if (!residues.empty()) {
            gemmi::Residue *res = residues[0];
            // Find parent chain for chain ID
            const char *chain_id = nullptr;
            for (auto &model : st.models) {
                for (auto &chain : model.chains) {
                    for (auto &r : chain.residues) {
                        if (&r == res) {
                            chain_id = chain.name.c_str();
                        }
                    }
                }
            }
            int seqnum = res->seqid.num.value;
            const char *resname = res->name.c_str();
            EXPECT_EQ(std::string(chain_id), "A");
            EXPECT_EQ(seqnum, 10);
            EXPECT_EQ(std::string(resname), "ILE");
        }
    }

    // Case 2: Multi-residue CID using "||" separator (OR logic)
    {
        std::string cid = "//A/10||//A/20";
        std::vector<gemmi::Residue *> residues = coot::cid_to_residues_gemmi(cid, st);
        EXPECT_EQ(residues.size(), 2u);
        for (size_t i = 0; i < residues.size(); ++i) {
            gemmi::Residue *res = residues[i];
            const char *chain_id = nullptr;
            for (auto &model : st.models) {
                for (auto &chain : model.chains) {
                    for (auto &r : chain.residues) {
                        if (&r == res) {
                            chain_id = chain.name.c_str();
                        }
                    }
                }
            }
            int seqnum = res->seqid.num.value;
            const char *resname = res->name.c_str();
            if (i == 0) {
                EXPECT_EQ(std::string(chain_id), "A");
                EXPECT_EQ(seqnum, 10);
                EXPECT_EQ(std::string(resname), "ILE");
            } else {
                EXPECT_EQ(std::string(chain_id), "A");
                EXPECT_EQ(seqnum, 20);
                EXPECT_EQ(std::string(resname), "LYS");
            }
        }
    }

    // Case 3: Non-existent residue (edge case — should return empty)
    {
        std::string cid = "//A/9999";
        std::vector<gemmi::Residue *> residues = coot::cid_to_residues_gemmi(cid, st);
        EXPECT_EQ(residues.size(), 0u);
    }

    // Case 4: Chain B residue (to verify multiple chains — chain B has only residue 1)
    {
        std::string cid = "//B/1";
        std::vector<gemmi::Residue *> residues = coot::cid_to_residues_gemmi(cid, st);
        EXPECT_EQ(residues.size(), 1u);
        if (!residues.empty()) {
            gemmi::Residue *res = residues[0];
            const char *chain_id = nullptr;
            for (auto &model : st.models) {
                for (auto &chain : model.chains) {
                    for (auto &r : chain.residues) {
                        if (&r == res) {
                            chain_id = chain.name.c_str();
                        }
                    }
                }
            }
            int seqnum = res->seqid.num.value;
            const char *resname = res->name.c_str();
            EXPECT_EQ(std::string(chain_id), "B");
            EXPECT_EQ(seqnum, 1);
            EXPECT_EQ(std::string(resname), "HOH");
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}