#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <fstream>
#include <string>
#include "function.hh"

bool file_exists(const std::string &path) {
    std::ifstream f(path);
    return f.good();
}

TEST(OracleTest, minimol_write_cif) {
    // Case 1: valid molecule from PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int atoms_loaded = static_cast<int>(gemmi::count_atom_sites(st));
        EXPECT_EQ(atoms_loaded, 2107);

        std::string cif_out = "/tmp/oracle_test_output.cif";
        int result = coot::minimol::molecule::write_cif_file_gemmi(st, cif_out);
        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists(cif_out));
    }

    // Case 2: empty molecule (edge case)
    {
        gemmi::Structure empty_st;
        int atoms = static_cast<int>(gemmi::count_atom_sites(empty_st));
        EXPECT_EQ(atoms, 0);

        std::string cif_out = "/tmp/oracle_test_empty.cif";
        int result = coot::minimol::molecule::write_cif_file_gemmi(empty_st, cif_out);
        EXPECT_EQ(result, 0);
        EXPECT_TRUE(file_exists(cif_out));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}