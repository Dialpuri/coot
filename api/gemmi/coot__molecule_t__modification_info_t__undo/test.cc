#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, UndoTest) {
  // Set up modification_info directly
  coot::modification_info_t::modification_info_t mod_info;
  mod_info.mol_name = "test-molecule";
  mod_info.is_mmcif_flag = false;
  
  // Load PDB using gemmi
  gemmi::Structure mol = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  
  // Record initial state
  EXPECT_EQ(mod_info.save_info.size(), 0u);
  EXPECT_EQ(mod_info.modification_index, 0);
  
  // Make initial backup
  std::string msg = mod_info.make_backup(mol, "initial");
  EXPECT_EQ(msg, "") << "make_backup should succeed with empty message";
  
  EXPECT_EQ(mod_info.save_info.size(), 1u);
  EXPECT_EQ(mod_info.modification_index, 1);
  
  // Call undo (using gemmi version)
  gemmi::Structure result = mod_info.undo_gemmi(mol);
  
  // Verify undo behavior: result should be non-empty (we have backup)
  EXPECT_FALSE(result.models.empty());
  
  // Result should have the same content as the backup (index 0)
  EXPECT_EQ(result.models.size(), mol.models.size());
  
  // Check that modification_index was updated (undid back to index 0)
  EXPECT_EQ(mod_info.modification_index, 0);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}