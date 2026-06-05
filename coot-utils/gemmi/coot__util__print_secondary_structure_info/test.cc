#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <iostream>
#include "function.hh"

TEST(GemmiSecStructureTest, print_secondary_structure_info) {
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-dev/coot/api/gideondoesntapprove.pdb");

  EXPECT_FALSE(st.models.empty());
  EXPECT_FALSE(st.models[0].chains.empty());

  std::cout << "Before: helices=" << st.helices.size()
            << " sheets=" << st.sheets.size() << std::endl;

  coot::util::print_secondary_structure_info_gemmi(st);

  // The PDB file has HELIX/SHEET/STRAND records, so we expect non-zero counts
  EXPECT_GT(st.helices.size(), 0) << "Expected some helices from PDB HELIX records";
  EXPECT_GT(st.sheets.size(), 0) << "Expected some sheets from PDB SHEET records";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}