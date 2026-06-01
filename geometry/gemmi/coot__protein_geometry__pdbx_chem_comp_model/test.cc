#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, pdbx_chem_comp_model_gemmi) {
  // Case 1: Create a Block with an "id" tag set to "LEU"
  {
    gemmi::cif::Block block("_chem_comp");
    block.set_pair("_chem_comp.id", "LEU");
    block.set_pair("_chem_comp.type", "L-LEUCINE");

    std::string result = coot::protein_geometry::pdbx_chem_comp_model_gemmi(block, 0);
    EXPECT_EQ(result, "LEU");
  }

  // Case 2: Create a Block with an "id" tag set to empty
  {
    gemmi::cif::Block block("_chem_comp");
    block.set_pair("_chem_comp.id", "");

    std::string result = coot::protein_geometry::pdbx_chem_comp_model_gemmi(block, 0);
    EXPECT_EQ(result, "");
  }

  // Case 3: No "id" tag present — should return empty string
  {
    gemmi::cif::Block block("_chem_comp");
    block.set_pair("_chem_comp.type", "SOME_TYPE");

    std::string result = coot::protein_geometry::pdbx_chem_comp_model_gemmi(block, 0);
    EXPECT_EQ(result, "");
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}