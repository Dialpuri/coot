#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_torsions_for_residues_in_chain_as_json) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid chain "A"
    {
        std::string chain_id = "A";
        std::string json = coot::molecule_t::get_torsions_for_residues_in_chain_as_json_gemmi(st, chain_id);
        EXPECT_EQ(json.length(), 57227u);
        EXPECT_FALSE(json.empty());
        EXPECT_NE(json.find("\"chain_id\":\"A\""), std::string::npos);
    }

    // Case 2: valid chain "B"
    {
        std::string chain_id = "B";
        std::string json = coot::molecule_t::get_torsions_for_residues_in_chain_as_json_gemmi(st, chain_id);
        EXPECT_EQ(json.length(), 2u);
        EXPECT_EQ(json, "[]");
    }

    // Case 3: non-existent chain "Z" (edge case)
    {
        std::string chain_id = "Z";
        std::string json = coot::molecule_t::get_torsions_for_residues_in_chain_as_json_gemmi(st, chain_id);
        EXPECT_EQ(json.length(), 2u);
        EXPECT_EQ(json, "[]");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}