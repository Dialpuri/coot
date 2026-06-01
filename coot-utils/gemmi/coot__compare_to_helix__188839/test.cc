#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, compare_to_helix) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    // Find chain A
    gemmi::Chain *chain_a = nullptr;
    for (gemmi::Chain &ch : model->chains) {
        if (ch.name == "A") {
            chain_a = &ch;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // --- Case 1: 4 consecutive residues in chain A (indices 50-53 → seqnums 60-63) ---
    {
        std::vector<gemmi::Residue *> res_4;
        for (int seq = 50; seq <= 53; seq++) {
            if (static_cast<size_t>(seq) < chain_a->residues.size()) {
                res_4.push_back(&chain_a->residues[seq]);
            }
        }
        EXPECT_EQ(res_4.size(), 4u);

        coot::helical_results_t result = coot::compare_to_helix_gemmi(res_4);
        EXPECT_FALSE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 11.3608, 1e-4);
    }

    // --- Case 2: 4 consecutive residues in chain A (indices 100-103 → seqnums 113-116) ---
    {
        std::vector<gemmi::Residue *> res_100;
        for (int seq = 100; seq <= 103; seq++) {
            if (static_cast<size_t>(seq) < chain_a->residues.size()) {
                res_100.push_back(&chain_a->residues[seq]);
            }
        }
        EXPECT_EQ(res_100.size(), 4u);

        coot::helical_results_t result = coot::compare_to_helix_gemmi(res_100);
        EXPECT_TRUE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 1.50364, 1e-4);
    }

    // --- Case 3: Wrong size (3 residues) — guard clause skips core logic ---
    {
        std::vector<gemmi::Residue *> res_3;
        for (int seq = 50; seq <= 52; seq++) {
            if (static_cast<size_t>(seq) < chain_a->residues.size()) {
                res_3.push_back(&chain_a->residues[seq]);
            }
        }
        EXPECT_EQ(res_3.size(), 3u);

        coot::helical_results_t result = coot::compare_to_helix_gemmi(res_3);
        EXPECT_FALSE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}