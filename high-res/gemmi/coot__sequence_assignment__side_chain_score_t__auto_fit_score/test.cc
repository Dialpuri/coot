#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, auto_fit_score) {
    coot::sequence_assignment::side_chain_score_t scorer;

    // Case 1: chain A, seq 10, score as ALA
    {
        std::string cid = "A";
        int seq = 10;
        auto idx = coot::sequence_assignment::side_chain_name_index::ALA;
        float score = scorer.auto_fit_score_gemmi(cid, seq, idx);
        EXPECT_NEAR(score, 0.0f, 1e-4);
    }

    // Case 2: chain A, seq 20, score as GLY
    {
        std::string cid = "A";
        int seq = 20;
        auto idx = coot::sequence_assignment::side_chain_name_index::GLY;
        float score = scorer.auto_fit_score_gemmi(cid, seq, idx);
        EXPECT_NEAR(score, 0.0f, 1e-4);
    }

    // Case 3: chain A, seq 0, score as TRP
    {
        std::string cid = "A";
        int seq = 0;
        auto idx = coot::sequence_assignment::side_chain_name_index::TRP;
        float score = scorer.auto_fit_score_gemmi(cid, seq, idx);
        EXPECT_NEAR(score, 0.0f, 1e-4);
    }

    // Case 4: empty chain id, seq 50, score as VAL
    {
        std::string cid = "";
        int seq = 50;
        auto idx = coot::sequence_assignment::side_chain_name_index::VAL;
        float score = scorer.auto_fit_score_gemmi(cid, seq, idx);
        EXPECT_NEAR(score, 0.0f, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}