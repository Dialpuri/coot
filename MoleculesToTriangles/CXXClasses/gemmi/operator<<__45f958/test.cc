#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <sstream>
#include "function.hh"

TEST(OracleTest, GetNumberOfAtoms) {
    // Case 1: Load PDB and check atom count
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int n_atoms = static_cast<int>(coot::molecule_t::get_number_of_atoms_gemmi(st));
        EXPECT_EQ(n_atoms, 2108);
    }

    // Case 2: Hydrogen PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        int n_atoms = static_cast<int>(coot::molecule_t::get_number_of_atoms_gemmi(st));
        EXPECT_EQ(n_atoms, 4281);
    }

    // Case 3: Empty manager (no file loaded)
    {
        gemmi::Structure st;  // empty structure

        int n_atoms = static_cast<int>(coot::molecule_t::get_number_of_atoms_gemmi(st));
        EXPECT_EQ(n_atoms, 0);
    }
}

TEST(OracleTest, OperatorLtgt) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    std::ostringstream oss;
    operator_ltlt_gemmi(oss, st);
    std::string output = oss.str();
    EXPECT_NE(output.find("Original name:"), std::string::npos);
    EXPECT_NE(output.find("nAtoms:"), std::string::npos);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}