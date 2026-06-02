#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, minimol_molecule_constructor) {
    // Load PDB via gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid structure pointer, udd_atom_index_to_user_data = false
    {
        coot::minimol::molecule m = coot::minimol::molecule_gemmi(&st, false);
        EXPECT_EQ(m.fragments.size(), 2u);
        EXPECT_TRUE(m.has_atoms());
        EXPECT_EQ(m.get_number_of_atoms(), 2107);
        EXPECT_EQ((int)(m.get_cell().size() > 0), 1);
        EXPECT_EQ((int)(!m.get_spacegroup().empty()), 1);
        EXPECT_EQ(m.mmdb_spacegroup, "P 21 21 21");
    }

    // Case 2: valid structure pointer, udd_atom_index_to_user_data = true
    {
        coot::minimol::molecule m = coot::minimol::molecule_gemmi(&st, true);
        EXPECT_EQ(m.fragments.size(), 2u);
        EXPECT_TRUE(m.has_atoms());
        EXPECT_EQ(m.get_number_of_atoms(), 2107);
    }

    // Case 3: NULL structure pointer (edge case — should produce empty molecule)
    {
        coot::minimol::molecule m = coot::minimol::molecule_gemmi(nullptr, false);
        EXPECT_EQ(m.fragments.size(), 0u);
        EXPECT_FALSE(m.has_atoms());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}