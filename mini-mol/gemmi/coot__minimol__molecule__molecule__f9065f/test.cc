#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__molecule__setup/gemmi/function.hh"

TEST(OracleTest, molecule) {
    // Case 1: Default constructor — empty molecule
    {
        coot::minimol::molecule mol = coot::minimol::molecule_gemmi();
        EXPECT_EQ(static_cast<short int>(mol.get_cell().size() == 6), 0);
        EXPECT_EQ(static_cast<short int>(!mol.get_spacegroup().empty()), 0);
        EXPECT_EQ(mol.fragments.size(), 0);
        EXPECT_TRUE(mol.is_empty());
        EXPECT_FALSE(mol.has_atoms());
    }

    // Case 2: After default construction, then set_cell and set_spacegroup
    {
        coot::minimol::molecule mol = coot::minimol::molecule_gemmi();
        float cell_vals[6] = {50.0f, 50.0f, 50.0f, 90.0f, 90.0f, 90.0f};
        mol.set_cell(cell_vals);
        mol.set_spacegroup("P 1");
        EXPECT_EQ(static_cast<short int>(mol.get_cell().size() == 6), 1);
        EXPECT_EQ(static_cast<short int>(!mol.get_spacegroup().empty()), 1);
        EXPECT_EQ(mol.get_cell().size(), 6u);
        EXPECT_EQ(mol.get_spacegroup(), "P 1");
    }

    // Case 3: Read from PDB file
    {
        coot::minimol::molecule mol = coot::minimol::molecule_gemmi();
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        short int result = coot::minimol::setup_gemmi(mol, &st);
        EXPECT_EQ(result, 0);
        EXPECT_EQ(mol.fragments.size(), 2u);
        EXPECT_TRUE(mol.has_atoms());
        EXPECT_EQ(mol.get_number_of_atoms(), 2107);
        EXPECT_EQ(static_cast<short int>(mol.get_cell().size() == 6), 1);
        EXPECT_EQ(static_cast<short int>(!mol.get_spacegroup().empty()), 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
