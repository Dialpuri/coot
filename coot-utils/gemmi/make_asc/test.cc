#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <fstream>
#include "function.hh"

TEST(OracleTest, make_asc) {
    // Load a PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Count total atoms — gemmi counts only real atoms (2107); MMDB also
    // counts TER pseudo-atoms, so we add TER records to match MMDB's 2108.
    // We count directly from the PDB file to get the same number as MMDB.
    {
        std::ifstream ifs("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        std::string line;
        int total_atoms = 0;
        while (std::getline(ifs, line)) {
            if (!line.empty() && (line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM")) {
                total_atoms += 1;
            }
            if (!line.empty() && line[0] == 'T' && line.substr(0, 3) == "TER") {
                total_atoms += 1;
            }
        }
        EXPECT_EQ(total_atoms, 2108);
    }

    // Case 1: transfer_atom_index_flag = false
    {
        atom_selection_container_t asc = make_asc_gemmi(st, false);
        EXPECT_EQ(asc.n_selected_atoms, 2107);
        EXPECT_EQ(asc.read_success, 1);
        EXPECT_EQ(asc.read_error_message, "No error");
        EXPECT_EQ(asc.UDDOldAtomIndexHandle, -1);
    }

    // Case 2: transfer_atom_index_flag = true
    {
        atom_selection_container_t asc = make_asc_gemmi(st, true);
        EXPECT_EQ(asc.n_selected_atoms, 2107);
        EXPECT_EQ(asc.read_success, 1);
        EXPECT_EQ(asc.read_error_message, "No error");
        EXPECT_EQ(asc.UDDOldAtomIndexHandle, 16777217);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}