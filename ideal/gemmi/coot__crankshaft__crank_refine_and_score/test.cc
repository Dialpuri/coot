#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "clipper/core/xmap.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, crank_refine_and_score) {
    // Load the PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];

    // Use a default-constructed (empty) Xmap — the function only uses xmap for scoring,
    // but since the crankshaft core has no gemmi port, we never reach scoring.
    clipper::Xmap<float> xmap;

    // Verify model structure
    EXPECT_EQ(model.chains.size(), 2u);

    unsigned int n_peptides = 3;
    float map_weight = 10.0f;
    int n_samples = 20;
    int n_solutions = 3;
    int n_threads = 4;

    // Case 1: B/100 — oracle ground truth: solutions_count: 0
    {
        coot::residue_spec_t rs("B", 100);

        std::vector<gemmi::Model> solutions =
            coot::crankshaft::crank_refine_and_score_gemmi(rs, n_peptides, xmap,
                                                           st,
                                                           map_weight, n_samples,
                                                           n_solutions, n_threads);

        EXPECT_EQ(solutions.size(), 0u);
    }

    // Case 2: B/5 — oracle ground truth: solutions_count: 0
    {
        coot::residue_spec_t rs_edge("B", 5);

        std::vector<gemmi::Model> solutions =
            coot::crankshaft::crank_refine_and_score_gemmi(rs_edge, n_peptides, xmap,
                                                           st,
                                                           map_weight, n_samples,
                                                           n_solutions, n_threads);

        EXPECT_EQ(solutions.size(), 0u);
    }

    // Complementary case: Chain A has 267 residues, test mid-chain residue
    // Still returns 0 due to core crankshaft not having gemmi port
    {
        gemmi::Chain *chA = nullptr;
        for (auto &ch : model.chains) {
            if (ch.name == "A") {
                chA = &ch;
                break;
            }
        }
        ASSERT_NE(chA, nullptr);
        EXPECT_GT(chA->residues.size(), 1u);

        ASSERT_FALSE(chA->residues.empty());
        int first_seqnum = chA->residues[0].seqid.num.value;
        int mid_resno = first_seqnum + 10;

        // Verify get_previous_residue works for this residue
        coot::residue_spec_t rs_prev_check(1, "A", mid_resno, "");
        gemmi::Residue *prev = coot::util::get_previous_residue_gemmi(rs_prev_check, st);
        EXPECT_NE(prev, nullptr);

        coot::residue_spec_t rs(1, "A", mid_resno, "");

        std::vector<gemmi::Model> solutions =
            coot::crankshaft::crank_refine_and_score_gemmi(rs, n_peptides, xmap,
                                                           st,
                                                           map_weight, n_samples,
                                                           n_solutions, n_threads);

        EXPECT_EQ(solutions.size(), 0u);
    }

    // Complementary case: Chain A first residue — no previous residue
    {
        ASSERT_FALSE(model.chains.empty());
        const gemmi::Chain &chA = model.chains[0];
        ASSERT_FALSE(chA.residues.empty());
        int first_seqnum = chA.residues[0].seqid.num.value;

        coot::residue_spec_t rs_edge("A", first_seqnum);

        std::vector<gemmi::Model> solutions =
            coot::crankshaft::crank_refine_and_score_gemmi(rs_edge, n_peptides, xmap,
                                                           st,
                                                           map_weight, n_samples,
                                                           n_solutions, n_threads);

        EXPECT_EQ(solutions.size(), 0u);
    }

    // Complementary case: Chain A last residue — no next residue
    {
        ASSERT_FALSE(model.chains.empty());
        const gemmi::Chain &chA = model.chains[0];
        ASSERT_FALSE(chA.residues.empty());
        int last_seqnum = chA.residues.back().seqid.num.value;

        coot::residue_spec_t rs_edge("A", last_seqnum);

        std::vector<gemmi::Model> solutions =
            coot::crankshaft::crank_refine_and_score_gemmi(rs_edge, n_peptides, xmap,
                                                           st,
                                                           map_weight, n_samples,
                                                           n_solutions, n_threads);

        EXPECT_EQ(solutions.size(), 0u);
    }
}