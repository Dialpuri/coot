#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, centre_of_residues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid residues (residues 10-12 from chain A)
    {
        std::vector<gemmi::Residue*> residues;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq >= 10 && seq <= 12) {
                            residues.push_back(&res);
                        }
                    }
                }
            }
        }
        EXPECT_EQ(residues.size(), 3u);
        auto result = coot::centre_of_residues_gemmi(residues);
        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second.x, 24.6157, 1e-4);
        EXPECT_NEAR(result.second.y, 13.9767, 1e-4);
        EXPECT_NEAR(result.second.z, 65.2359, 1e-4);
    }

    // Case 2: empty vector (edge case — no atoms)
    {
        std::vector<gemmi::Residue*> residues;
        EXPECT_EQ(residues.size(), 0u);
        auto result = coot::centre_of_residues_gemmi(residues);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second.x, 0.0, 1e-4);
        EXPECT_NEAR(result.second.y, 0.0, 1e-4);
        EXPECT_NEAR(result.second.z, 0.0, 1e-4);
    }

    // Case 3: invalid residue CID (should result in empty vector)
    {
        std::vector<gemmi::Residue*> residues;
        // Search for a non-existent residue (seq 99999 in chain A)
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name == "A") {
                    for (auto& res : chain.residues) {
                        int seq = res.seqid.num.value;
                        if (seq == 99999) {
                            residues.push_back(&res);
                        }
                    }
                }
            }
        }
        EXPECT_EQ(residues.size(), 0u);
        auto result = coot::centre_of_residues_gemmi(residues);
        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second.x, 0.0, 1e-4);
        EXPECT_NEAR(result.second.y, 0.0, 1e-4);
        EXPECT_NEAR(result.second.z, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}