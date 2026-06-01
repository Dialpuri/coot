#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, side_chain_score) {
    coot::sequence_assignment::side_chain_score_t scorer;

    // cache_standard_residues populates the standard_residues vector
    short int rc = scorer.cache_standard_residues();
    EXPECT_EQ(rc, 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}