#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, have_intermediate_residue_by_seqnum) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Find chain A
    gemmi::Chain* chainA = nullptr;
    gemmi::Chain* chainB = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") chainA = &chain;
        if (chain.name == "B") chainB = &chain;
    }
    ASSERT_NE(chainA, nullptr);

    // Build residues vector with all residues from chain A
    std::vector<std::pair<bool, gemmi::CRA>> residues_vec;
    for (gemmi::Residue& res : chainA->residues) {
        gemmi::CRA cra{chainA, &res, nullptr};
        residues_vec.push_back(std::make_pair(false, cra));
    }

    // Also add a residue from chain B for cross-chain test
    if (chainB) {
        gemmi::Residue& resB = chainB->residues[0];
        gemmi::CRA craB{chainB, &resB, nullptr};
        residues_vec.push_back(std::make_pair(false, craB));
    }

    // Case 1: Adjacent residues (seqnum diff = 1) -> should return false
    {
        ASSERT_GE(chainA->residues.size(), 2u);
        gemmi::Residue& res1 = chainA->residues[0];   // seqnum 1
        gemmi::Residue& res2 = chainA->residues[1];   // seqnum 2
        gemmi::CRA cra1{chainA, &res1, nullptr};
        gemmi::CRA cra2{chainA, &res2, nullptr};
        bool result = coot::have_intermediate_residue_by_seqnum_gemmi(cra1, cra2, residues_vec);
        EXPECT_EQ(result, false);
    }

    // Case 2: Non-adjacent residues with intermediate in vec (seqnum diff > 1) -> should return true
    {
        ASSERT_GE(chainA->residues.size(), 3u);
        gemmi::Residue& first = chainA->residues[0];   // seqnum 1
        gemmi::Residue& third = chainA->residues[2];   // seqnum 3
        gemmi::CRA cra1{chainA, &first, nullptr};
        gemmi::CRA cra3{chainA, &third, nullptr};
        bool result = coot::have_intermediate_residue_by_seqnum_gemmi(cra1, cra3, residues_vec);
        EXPECT_EQ(result, true);
    }

    // Case 3: Residues from different chains -> should return false
    {
        if (chainB && !chainB->residues.empty()) {
            gemmi::Residue& res1 = chainA->residues[0];
            gemmi::Residue& resB = chainB->residues[0];
            gemmi::CRA cra1{chainA, &res1, nullptr};
            gemmi::CRA craB{chainB, &resB, nullptr};
            bool result = coot::have_intermediate_residue_by_seqnum_gemmi(cra1, craB, residues_vec);
            EXPECT_EQ(result, false);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}