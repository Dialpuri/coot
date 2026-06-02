#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, minimol_molecule_setup) {
    // Case 1: valid molecule from a real PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        EXPECT_EQ(st.models.size(), 1u);

        coot::minimol::molecule mm;

        // State BEFORE setup
        EXPECT_EQ(mm.fragments.size(), 0);
        EXPECT_EQ(static_cast<int>(mm.get_cell().size()), 0);
        EXPECT_TRUE(mm.get_spacegroup().empty());

        short istat = coot::minimol::setup_gemmi(mm, &st);

        EXPECT_EQ(static_cast<int>(istat), 0);
        EXPECT_EQ(mm.fragments.size(), 2);
        EXPECT_EQ(static_cast<int>(mm.get_cell().size()), 6);
        EXPECT_FALSE(mm.get_spacegroup().empty());
        EXPECT_EQ(mm.get_spacegroup(), "P 21 21 21");
        EXPECT_EQ(mm.mmdb_cell.size(), 6);
        EXPECT_EQ(mm.get_number_of_atoms(), 2107);
    }

    // Case 2: NULL manager — triggers the guard clause
    {
        coot::minimol::molecule mm2;
        EXPECT_EQ(mm2.fragments.size(), 0);

        short istat2 = coot::minimol::setup_gemmi(mm2, nullptr);

        EXPECT_EQ(static_cast<int>(istat2), 1);
        EXPECT_EQ(mm2.fragments.size(), 0);
        EXPECT_EQ(static_cast<int>(mm2.get_cell().size()), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}