#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a gemmi::Residue by chain name and seqnum, return CRA
static gemmi::CRA find_residue(gemmi::Structure& st,
                                const std::string& chain_name,
                                int seqnum) {
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
    }
    return gemmi::CRA{nullptr, nullptr, nullptr};
}

TEST(OracleTest, bonded_pair_reorder_as_needed) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::CRA cra10 = find_residue(st, "A", 10);  // ILE 10
    gemmi::CRA cra5  = find_residue(st, "A", 5);   // GLN 5
    gemmi::CRA cra50 = find_residue(st, "A", 50);  // ARG 50
    ASSERT_NE(cra10.residue, nullptr);
    ASSERT_NE(cra5.residue, nullptr);
    ASSERT_NE(cra50.residue, nullptr);

    // Case 1: res_2 seqnum < res_1 seqnum, same chain, both amino acids -> swap
    {
        coot::bonded_pair_t bp(cra10, cra5, false, false, "peptide");
        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 10);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 5);
        EXPECT_EQ(bp.is_fixed_first, false);
        EXPECT_EQ(bp.is_fixed_second, false);

        bp.reorder_as_needed_gemmi();

        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 5);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 10);
        EXPECT_EQ(bp.is_fixed_first, false);
        EXPECT_EQ(bp.is_fixed_second, false);
    }

    // Case 2: already correctly ordered -> no swap
    {
        coot::bonded_pair_t bp(cra5, cra10, false, false, "peptide");
        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 5);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 10);
        EXPECT_EQ(bp.is_fixed_first, false);
        EXPECT_EQ(bp.is_fixed_second, false);

        bp.reorder_as_needed_gemmi();

        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 5);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 10);
        EXPECT_EQ(bp.is_fixed_first, false);
        EXPECT_EQ(bp.is_fixed_second, false);
    }

    // Case 3: is_fixed_first=true, is_fixed_second=false with swap needed
    // Both fixed flags swap alongside residue pointers
    {
        coot::bonded_pair_t bp(cra50, cra5, true, false, "contact");
        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 50);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 5);
        EXPECT_EQ(bp.is_fixed_first, true);
        EXPECT_EQ(bp.is_fixed_second, false);

        bp.reorder_as_needed_gemmi();

        EXPECT_STREQ(bp.res_1.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_1.residue->seqid.num.value, 5);
        EXPECT_STREQ(bp.res_2.chain->name.c_str(), "A");
        EXPECT_EQ(bp.res_2.residue->seqid.num.value, 50);
        EXPECT_EQ(bp.is_fixed_first, false);
        EXPECT_EQ(bp.is_fixed_second, true);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}