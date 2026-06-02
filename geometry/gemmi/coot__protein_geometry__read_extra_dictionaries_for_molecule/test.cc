#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, read_extra_dictionaries_for_molecule) {
    // Load PDB into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Build protein geometry with standard dictionaries
    coot::protein_geometry geom;
    geom.init_standard();

    int imol = 0;

    // Case 1: Normal call with valid molecule (standard protein - all should have dicts)
    {
        int read_number = 0;
        bool success = coot::read_extra_dictionaries_for_molecule_gemmi(geom, &st, imol, &read_number);
        EXPECT_TRUE(success);
        EXPECT_EQ(read_number, 0);
    }

    // Case 2: Null mol pointer — guarded path
    {
        int read_number = 0;
        bool success = coot::read_extra_dictionaries_for_molecule_gemmi(geom, nullptr, imol, &read_number);
        EXPECT_FALSE(success);
        EXPECT_EQ(read_number, 0);
    }

    // Case 3: Different imol value (should still find no missing dicts for standard residues)
    {
        int read_number = 42;
        bool success = coot::read_extra_dictionaries_for_molecule_gemmi(geom, &st, 5, &read_number);
        EXPECT_TRUE(success);
        EXPECT_EQ(read_number, 42);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}