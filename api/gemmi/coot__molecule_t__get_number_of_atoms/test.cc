#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, molecules_container_t_read_pdb) {
    const std::string pdb_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

    // Edge case 1: read a valid PDB
    gemmi::Structure st = gemmi::read_pdb_file(pdb_file);
    unsigned int n1 = coot::molecule_t::get_number_of_atoms_gemmi(st);
    EXPECT_EQ(n1, 2107u);

    // Edge case 2: try reading a non-existent file
    // gemmi::read_pdb_file will throw or return an empty structure
    try {
        gemmi::Structure st2 = gemmi::read_pdb_file("/nonexistent/file.pdb");
        // If it doesn't throw, check if it's empty
        unsigned int n2 = coot::molecule_t::get_number_of_atoms_gemmi(st2);
        EXPECT_EQ(n2, 0u);
    } catch (...) {
        // Expected - file doesn't exist
    }

    // Edge case 3: call on molecule 0 (if valid)
    unsigned int n3 = coot::molecule_t::get_number_of_atoms_gemmi(st);
    EXPECT_EQ(n3, 2107u);

    // Edge case 4: check number of molecules loaded
    // In gemmi, we just have one structure
    EXPECT_EQ(st.models.size(), 1u);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}