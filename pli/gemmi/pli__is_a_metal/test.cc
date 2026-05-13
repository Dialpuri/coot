#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, is_a_metal) {
  // Load the PDB file
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  ASSERT_FALSE(st.models.empty());

  // Case 1: Test with ILE (non-metal residue)
  {
    bool found_ile = false;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        for (auto& res : chain.residues)
          if (res.name == "ILE" && res.seqid.num.value == 10) {
            found_ile = true;
            bool result = pli::is_a_metal_gemmi(res);
            EXPECT_EQ(result, false);
          }
    EXPECT_TRUE(found_ile) << "Residue ILE at seqnum 10 not found";
  }

  // Case 2: Test with LYS (non-metal residue)
  {
    bool found_lys = false;
    for (auto& model : st.models)
      for (auto& chain : model.chains)
        for (auto& res : chain.residues)
          if (res.name == "LYS" && res.seqid.num.value == 20) {
            found_lys = true;
            bool result = pli::is_a_metal_gemmi(res);
            EXPECT_EQ(result, false);
          }
    EXPECT_TRUE(found_lys) << "Residue LYS at seqnum 20 not found";
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}