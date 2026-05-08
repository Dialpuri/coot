#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chain_info_tConstructor) {
  // Load the PDB file
  const char *pdb_file = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
  gemmi::Structure st = gemmi::read_pdb_file(pdb_file);

  // Get model 1 (gemmi is 0-indexed, so model 1 → index 0)
  gemmi::Model &model = st.models[0];
  ASSERT_FALSE(model.chains.empty());

  // Get chain A (first chain)
  gemmi::Chain &chain = model.chains[0];

  // Find min and max residue numbers
  int resno_min = 999999;
  int resno_max = -999999;
  for (const auto &res : chain.residues) {
    int seq_num = res.seqid.num.value;
    if (seq_num < resno_min) resno_min = seq_num;
    if (seq_num > resno_max) resno_max = seq_num;
  }

  // Test the constructor with actual values from PDB
  auto info = coot::molecule_t::get_molecular_representation_mesh::chain_info_t_gemmi(&chain, resno_min, resno_max);

  // Verify the constructor worked
  EXPECT_EQ(info.chain_p, &chain);
  EXPECT_EQ(info.resno_min, resno_min);
  EXPECT_EQ(info.resno_max, resno_max);

  // Edge case: identical min/max
  {
    auto info2 = coot::molecule_t::get_molecular_representation_mesh::chain_info_t_gemmi(&chain, 10, 10);
    EXPECT_EQ(info2.chain_p, &chain);
    EXPECT_EQ(info2.resno_min, 10);
    EXPECT_EQ(info2.resno_max, 10);
  }

  // Edge case: negative min
  {
    auto info3 = coot::molecule_t::get_molecular_representation_mesh::chain_info_t_gemmi(&chain, -5, 20);
    EXPECT_EQ(info3.chain_p, &chain);
    EXPECT_EQ(info3.resno_min, -5);
    EXPECT_EQ(info3.resno_max, 20);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}