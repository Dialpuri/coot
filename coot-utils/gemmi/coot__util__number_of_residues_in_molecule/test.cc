#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, NumberOfResiduesInMolecule) {
    // case 1: valid structure (should count residues across all models/chains)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int n_res = coot::util::number_of_residues_in_molecule_gemmi(st);
        EXPECT_EQ(n_res, 268);
    }

    // case 2: empty structure (guard clause)
    {
        gemmi::Structure st;
        int n_res = coot::util::number_of_residues_in_molecule_gemmi(st);
        EXPECT_EQ(n_res, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}