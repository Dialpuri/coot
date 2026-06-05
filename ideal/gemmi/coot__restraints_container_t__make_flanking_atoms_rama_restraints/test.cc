#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_flanking_atoms_rama_restraints) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    // Case 1: Single middle residue with flanking neighbors (residue 10 in chain A)
    {
        const std::string chain_id_save = "A";
        int istart_res = 10;
        int iend_res = 10;
        int istart_minus_flag = 1;
        int iend_plus_flag = 1;

        auto triples = coot::find_rama_flanking_triples_gemmi(
            model, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag);

        EXPECT_EQ(1U, triples.size());
        if (!triples.empty()) {
            EXPECT_EQ(9, triples[0].first->seqid.num.value);
            EXPECT_EQ(10, triples[0].middle->seqid.num.value);
            EXPECT_EQ(11, triples[0].third->seqid.num.value);
        }

        EXPECT_EQ("A", chain_id_save);
        EXPECT_EQ(10, istart_res);
        EXPECT_EQ(10, iend_res);
        EXPECT_EQ(1, istart_minus_flag);
        EXPECT_EQ(1, iend_plus_flag);

        coot::protein_geometry pg;
        int result = coot::make_flanking_atoms_rama_restraints_gemmi(
            st, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag, pg);
        EXPECT_EQ(0, result);
    }

    // Case 2: Two moving residues (residues 10-11 in chain A)
    {
        const std::string chain_id_save = "A";
        int istart_res = 10;
        int iend_res = 11;
        int istart_minus_flag = 1;
        int iend_plus_flag = 1;

        auto triples = coot::find_rama_flanking_triples_gemmi(
            model, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag);

        EXPECT_EQ(2U, triples.size());
        if (triples.size() >= 2) {
            EXPECT_EQ(9, triples[0].first->seqid.num.value);
            EXPECT_EQ(10, triples[0].middle->seqid.num.value);
            EXPECT_EQ(11, triples[0].third->seqid.num.value);
            EXPECT_EQ(10, triples[1].first->seqid.num.value);
            EXPECT_EQ(11, triples[1].middle->seqid.num.value);
            EXPECT_EQ(12, triples[1].third->seqid.num.value);
        }

        EXPECT_EQ("A", chain_id_save);
        EXPECT_EQ(10, istart_res);
        EXPECT_EQ(11, iend_res);
        EXPECT_EQ(1, istart_minus_flag);
        EXPECT_EQ(1, iend_plus_flag);

        coot::protein_geometry pg;
        int result = coot::make_flanking_atoms_rama_restraints_gemmi(
            st, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag, pg);
        EXPECT_EQ(0, result);
    }

    // Case 3: Guard path — istart_minus_flag not set (should add nothing)
    {
        const std::string chain_id_save = "A";
        int istart_res = 10;
        int iend_res = 10;
        int istart_minus_flag = 0;
        int iend_plus_flag = 0;

        auto triples = coot::find_rama_flanking_triples_gemmi(
            model, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag);

        EXPECT_EQ(0U, triples.size());

        EXPECT_EQ("A", chain_id_save);
        EXPECT_EQ(10, istart_res);
        EXPECT_EQ(10, iend_res);
        EXPECT_EQ(0, istart_minus_flag);
        EXPECT_EQ(0, iend_plus_flag);

        coot::protein_geometry pg;
        int result = coot::make_flanking_atoms_rama_restraints_gemmi(
            st, chain_id_save, istart_res, iend_res,
            istart_minus_flag, iend_plus_flag, pg);
        EXPECT_EQ(0, result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}