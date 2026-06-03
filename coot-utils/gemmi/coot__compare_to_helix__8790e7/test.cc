#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, compare_to_helix) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain *chain = &model->chains[0];
    ASSERT_NE(chain, nullptr);

    // Gather consecutive residues
    std::vector<gemmi::Residue *> residues;
    for (auto &res : chain->residues) {
        residues.push_back(&res);
    }

    // Case 1: 4 consecutive residues (valid alpha helix test)
    {
        std::vector<gemmi::Residue *> helical_residues;
        for (int i = 0; i < 4; i++) {
            helical_residues.push_back(residues[i]);
        }

        coot::helical_results_t result = coot::compare_to_helix_gemmi(helical_residues);

        EXPECT_FALSE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 10.3456, 1e-4);
    }

    // Case 2: only 3 residues — should return empty/default results
    {
        std::vector<gemmi::Residue *> short_residues;
        for (int i = 0; i < 3; i++) {
            short_residues.push_back(residues[i]);
        }

        coot::helical_results_t result = coot::compare_to_helix_gemmi(short_residues);

        EXPECT_FALSE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 0.0, 1e-4);
    }

    // Case 3: empty vector — edge case
    {
        std::vector<gemmi::Residue *> empty_residues;
        coot::helical_results_t result = coot::compare_to_helix_gemmi(empty_residues);

        EXPECT_FALSE(result.is_alpha_helix_like);
        EXPECT_NEAR(result.sum_delta, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}