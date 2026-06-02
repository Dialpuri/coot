#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cache_standard_residues) {
    // Case 1: default-constructed scorer
    {
        coot::sequence_assignment::side_chain_score_t scorer;

        EXPECT_EQ(scorer.standard_residues.size(), 0u);

        short int result = scorer.cache_standard_residues_gemmi(nullptr);

        EXPECT_EQ(result, 1);
        EXPECT_EQ(scorer.standard_residues.size(), 0u);
    }

    // Case 2: second scorer instance (function does NOT overwrite existing data)
    {
        coot::sequence_assignment::side_chain_score_t scorer2;

        EXPECT_EQ(scorer2.standard_residues.size(), 0u);

        short int result2 = scorer2.cache_standard_residues_gemmi(nullptr);

        EXPECT_EQ(result2, 1);
        EXPECT_EQ(scorer2.standard_residues.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}