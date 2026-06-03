#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, chem_comp) {
  int imol_enc = -999999; // coot::protein_geometry::IMOL_ENC_ANY

  // Case 1: loop with one LZA entry
  {
    gemmi::cif::Loop loop;
    loop.tags = {"_chem_comp.id", "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
                 "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh",
                 "_chem_comp.release_status", "_chem_comp.desc_level"};
    loop.values = {"LZA", "LZA", "LZA_ligand", "L-peptide", "7", "2", "REL", "."};

    EXPECT_EQ(imol_enc, -999999);
    EXPECT_EQ(loop.length(), 1);

    std::string result = coot::protein_geometry::chem_comp_gemmi(loop, imol_enc);

    EXPECT_EQ(result, "LZA");
    EXPECT_EQ(loop.length(), 1);
  }

  // Case 2: empty loop — should return empty string
  {
    gemmi::cif::Loop empty_loop;
    empty_loop.tags = {"_chem_comp.id", "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
                       "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh",
                       "_chem_comp.release_status", "_chem_comp.desc_level"};
    // no values = 0 rows

    EXPECT_EQ(imol_enc, -999999);
    EXPECT_EQ(empty_loop.length(), 0);

    std::string result2 = coot::protein_geometry::chem_comp_gemmi(empty_loop, imol_enc);

    EXPECT_EQ(result2, "");
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}