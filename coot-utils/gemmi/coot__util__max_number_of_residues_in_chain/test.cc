#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, max_number_of_residues_in_chain) {
    // case 1: valid structure with 2 chains (A: 267 residues, B: 1 residue)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int max_residues = coot::util::max_number_of_residues_in_chain_gemmi(st);
        EXPECT_EQ(max_residues, 267);
    }

    // case 2: empty structure — should return -1
    {
        gemmi::Structure st;
        int max_residues = coot::util::max_number_of_residues_in_chain_gemmi(st);
        EXPECT_EQ(max_residues, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}