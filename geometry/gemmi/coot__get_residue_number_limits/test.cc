#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include <string>
#include "function.hh"

TEST(OracleTest, get_residue_number_limits) {
    // case 1: valid PDB with real chains
    {
        std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
        std::map<std::string, std::pair<int, int> > limits = coot::get_residue_number_limits_gemmi(st);
        EXPECT_EQ(limits.size(), 2u);

        std::map<std::string, std::pair<int, int> >::iterator it;

        it = limits.find("A");
        ASSERT_NE(it, limits.end());
        EXPECT_EQ(it->second.first, 1);
        EXPECT_EQ(it->second.second, 298);

        it = limits.find("B");
        ASSERT_NE(it, limits.end());
        EXPECT_EQ(it->second.first, 1);
        EXPECT_EQ(it->second.second, 1);
    }

    // case 2: empty structure (no models) — should return empty map
    {
        gemmi::Structure empty_st;
        std::map<std::string, std::pair<int, int> > empty_limits = coot::get_residue_number_limits_gemmi(empty_st);
        EXPECT_EQ(empty_limits.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}