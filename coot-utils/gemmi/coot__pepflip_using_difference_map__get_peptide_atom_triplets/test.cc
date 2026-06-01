#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, pepflip_using_difference_map_get_peptide_atom_triplets_gemmi) {
    // Case 1: standard call — should find peptide triplets in chain A
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1u);
        EXPECT_EQ(st.models[0].chains.size(), 2u);

        std::vector<flip_atom_triplet_t_gemmi> triplets =
            coot::pepflip_using_difference_map::get_peptide_atom_triplets_gemmi(st);
        EXPECT_EQ(triplets.size(), 263u);

        // Check first 3 triplets
        // Note: gemmi atom names are unpadded ("CA", "O") vs MMDB padded (" CA ", " O  ")
        {
            flip_atom_triplet_t_gemmi& t = triplets[0];
            EXPECT_STREQ(t.CA_this->name.c_str(), "CA");
            EXPECT_STREQ(t.O_this->name.c_str(), "O");
            EXPECT_STREQ(t.CA_next->name.c_str(), "CA");
            EXPECT_EQ(t.res_this->seqid.num.value, 1);
            EXPECT_EQ(t.res_next->seqid.num.value, 2);
        }
        {
            flip_atom_triplet_t_gemmi& t = triplets[1];
            EXPECT_STREQ(t.CA_this->name.c_str(), "CA");
            EXPECT_STREQ(t.O_this->name.c_str(), "O");
            EXPECT_STREQ(t.CA_next->name.c_str(), "CA");
            EXPECT_EQ(t.res_this->seqid.num.value, 2);
            EXPECT_EQ(t.res_next->seqid.num.value, 3);
        }
        {
            flip_atom_triplet_t_gemmi& t = triplets[2];
            EXPECT_STREQ(t.CA_this->name.c_str(), "CA");
            EXPECT_STREQ(t.O_this->name.c_str(), "O");
            EXPECT_STREQ(t.CA_next->name.c_str(), "CA");
            EXPECT_EQ(t.res_this->seqid.num.value, 3);
            EXPECT_EQ(t.res_next->seqid.num.value, 4);
        }
    }

    // Case 2: empty molecule — no model
    {
        gemmi::Structure empty_st;
        std::vector<flip_atom_triplet_t_gemmi> triplets =
            coot::pepflip_using_difference_map::get_peptide_atom_triplets_gemmi(empty_st);

        EXPECT_EQ(empty_st.models.size(), 0u);
        EXPECT_EQ(triplets.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}