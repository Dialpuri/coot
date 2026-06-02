#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Case 1: Load PDB and compute CO orientations
TEST(CO_orientations_gemmi_Test, CO_orientations_example_pdb) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Original test used mol->GetModel(1) which is 1-indexed;
    // gemmi models are 0-indexed, so we use st.models[0].
    if (st.models.empty()) {
        EXPECT_TRUE(true); // no models, nothing to test
        return;
    }

    auto scores = coot::util::CO_orientations_gemmi(st.models[0]);

    // Oracle output: score_count: 0
    EXPECT_EQ(scores.size(), 0);
}

// Case 2: empty mol (no models) - should produce 0 scores
TEST(CO_orientations_gemmi_Test, CO_orientations_empty) {
    gemmi::Model empty_model;

    auto empty_scores = coot::util::CO_orientations_gemmi(empty_model);

    // Oracle output: empty_score_count: 0
    EXPECT_EQ(empty_scores.size(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}