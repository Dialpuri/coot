#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, select_residues_modes) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Case 1: SINGLE mode — returns exact match only
    {
        coot::residue_spec_t spec("A", 10, "");
        std::vector<gemmi::Residue*> rv = coot::select_residues_gemmi(spec, "SINGLE", st);
        EXPECT_EQ(rv.size(), 1u);
        gemmi::Residue *res = rv[0];
        EXPECT_NE(res, nullptr);
        // Check chain by searching in the structure (gemmi Residue has no parent pointer)
        bool found_correct_chain = false;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& r : chain.residues) {
                    if (&r == res && chain.name == "A") {
                        found_correct_chain = true;
                    }
                }
            }
        }
        EXPECT_TRUE(found_correct_chain);
        EXPECT_EQ(res->seqid.num.value, 10);
    }

    // Case 2: TRIPLE mode — returns previous, current, next residues
    {
        coot::residue_spec_t spec("A", 10, "");
        std::vector<gemmi::Residue*> rv = coot::select_residues_gemmi(spec, "TRIPLE", st);
        EXPECT_EQ(rv.size(), 3u);
        EXPECT_EQ(rv[0]->seqid.num.value, 9);
        EXPECT_EQ(rv[1]->seqid.num.value, 10);
        EXPECT_EQ(rv[2]->seqid.num.value, 11);
    }

    // Case 3: LITERAL mode — converts to SINGLE internally, same result
    {
        coot::residue_spec_t spec("A", 10, "");
        std::vector<gemmi::Residue*> rv = coot::select_residues_gemmi(spec, "LITERAL", st);
        EXPECT_EQ(rv.size(), 1u);
    }

    // Case 4: CHAIN mode — returns all residues in the same chain
    {
        coot::residue_spec_t spec("A", 10, "");
        std::vector<gemmi::Residue*> rv = coot::select_residues_gemmi(spec, "CHAIN", st);
        EXPECT_EQ(rv.size(), 267u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}