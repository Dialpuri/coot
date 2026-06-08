#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <utility>

TEST(OracleTest, get_middle_pos) {
    // Load PDB into gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: call get_middle_pos with the loaded structure
    {
        std::pair<clipper::Coord_orth, gemmi::Structure*> result =
            coot::high_res::get_middle_pos_gemmi(st);

        EXPECT_EQ(result.first.format(), "xyz = (      34.9,     4.631,     45.93)");
        EXPECT_NE(result.second, nullptr);
    }

    // Case 2: empty structure — edge case with no atoms
    {
        gemmi::Structure empty_st;
        std::pair<clipper::Coord_orth, gemmi::Structure*> result =
            coot::high_res::get_middle_pos_gemmi(empty_st);

        EXPECT_NE(result.second, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
