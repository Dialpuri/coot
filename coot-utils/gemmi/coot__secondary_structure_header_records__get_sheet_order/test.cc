#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_sheet_order) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];

    // Find chain A
    gemmi::Chain *chain_a = nullptr;
    for (gemmi::Chain &ch : model.chains) {
        if (ch.name == "A") {
            chain_a = &ch;
            break;
        }
    }
    EXPECT_NE(chain_a, nullptr);

    // Verify inputs match oracle
    EXPECT_EQ(model.num, 1);
    EXPECT_EQ(chain_a->name, "A");

    // Collect ALL residues from chain A
    std::vector<gemmi::Residue*> all_residues;
    for (gemmi::Residue &res : chain_a->residues) {
        all_residues.push_back(&res);
    }

    // Strand 1: residues 20-25
    std::vector<gemmi::Residue*> strand1;
    for (gemmi::Residue *r : all_residues) {
        int seq = r->seqid.num.value;
        if (seq >= 20 && seq <= 25) {
            strand1.push_back(r);
        }
    }

    // Strand 2: residues 30-35
    std::vector<gemmi::Residue*> strand2;
    for (gemmi::Residue *r : all_residues) {
        int seq = r->seqid.num.value;
        if (seq >= 30 && seq <= 35) {
            strand2.push_back(r);
        }
    }

    // Strand 3: residues 40-45
    std::vector<gemmi::Residue*> strand3;
    for (gemmi::Residue *r : all_residues) {
        int seq = r->seqid.num.value;
        if (seq >= 40 && seq <= 45) {
            strand3.push_back(r);
        }
    }

    std::vector<std::vector<gemmi::Residue*>> strands_with_residues;
    if (!strand1.empty()) strands_with_residues.push_back(strand1);
    if (!strand2.empty()) strands_with_residues.push_back(strand2);
    if (!strand3.empty()) strands_with_residues.push_back(strand3);

    // Verify input: n_strands == 2 (strand3 is empty)
    EXPECT_EQ(strands_with_residues.size(), 2u);

    // Verify strand details
    if (!strand1.empty()) {
        EXPECT_EQ(strand1[0]->seqid.num.value, 20);
        EXPECT_EQ(strand1.back()->seqid.num.value, 25);
        EXPECT_EQ(strand1.size(), 6u);
    }
    if (!strand2.empty()) {
        EXPECT_EQ(strand2[0]->seqid.num.value, 30);
        EXPECT_EQ(strand2.back()->seqid.num.value, 35);
        EXPECT_EQ(strand2.size(), 6u);
    }

    // Case 1: Normal strands
    {
        auto result = get_sheet_order_gemmi(st, model, strands_with_residues);

        EXPECT_EQ(result.size(), 1u);
        if (!result.empty()) {
            EXPECT_EQ(result[0].size(), 2u);

            if (result[0].size() >= 1) {
                EXPECT_EQ(result[0][0].strand_idx, 0);
                EXPECT_EQ(static_cast<int>(result[0][0].sense), 0); // FIRST
            }
            if (result[0].size() >= 2) {
                EXPECT_EQ(result[0][1].strand_idx, 1);
                EXPECT_EQ(static_cast<int>(result[0][1].sense), 2); // ANTI_PARALLEL
            }
        }
    }

    // Case 2: Empty strands input
    {
        std::vector<std::vector<gemmi::Residue*>> empty_strands;
        auto result_empty = get_sheet_order_gemmi(st, model, empty_strands);
        EXPECT_EQ(result_empty.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}