#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_single_letter_codes_for_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model &model = st.models[0];

    // Case 1: valid protein chain "A" (267 residues)
    {
        std::string chain_id = "A";
        auto v = coot::get_single_letter_codes_for_chain_gemmi(model, chain_id);
        EXPECT_EQ(v.size(), 267u);
        EXPECT_FALSE(v.empty());
        EXPECT_EQ(v[0].first.chain_id, "A");
        EXPECT_EQ(v[0].first.res_no, 1);
        EXPECT_EQ(v[0].second, "M");
        EXPECT_EQ(v.back().first.chain_id, "A");
        EXPECT_EQ(v.back().first.res_no, 298);
        EXPECT_EQ(v.back().second, "L");
    }

    // Case 2: chain "B" — single HOH (water) residue
    {
        std::string chain_id = "B";
        auto v = coot::get_single_letter_codes_for_chain_gemmi(model, chain_id);
        EXPECT_EQ(v.size(), 1u);
        EXPECT_FALSE(v.empty());
        EXPECT_EQ(v[0].first.chain_id, "B");
        EXPECT_EQ(v[0].first.res_no, 1);
        EXPECT_EQ(v[0].second, "~");
    }

    // Case 3: non-existent chain "Z" — should return empty
    {
        std::string chain_id = "Z";
        auto v = coot::get_single_letter_codes_for_chain_gemmi(model, chain_id);
        EXPECT_EQ(v.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}