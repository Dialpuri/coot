#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test the redo_gemmi function
TEST(OracleTest, ModificationInfoRedo) {
   // Load a PDB file to create a valid gemmi::Structure
   const char* pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
   
   // Create save_info_t objects with gemmi structures
   std::vector<coot::molecule_t::modification_info_t::save_info_t> save_info;
   
   // Create save_info entries by reading PDB
   save_info.emplace_back(pdb_path, "initial state");
   save_info[0].structure = std::make_unique<gemmi::Structure>(gemmi::read_pdb_file(pdb_path));
   
   save_info.emplace_back(pdb_path, "after state 1");
   save_info[1].structure = std::make_unique<gemmi::Structure>(gemmi::read_pdb_file(pdb_path));
   
   save_info.emplace_back(pdb_path, "after state 2");
   save_info[2].structure = std::make_unique<gemmi::Structure>(gemmi::read_pdb_file(pdb_path));
   
   // Set modification_index to a position where redo() can find something
   // After 3 backups, modification_index == 3 (save_info.size())
   // We set it to 1 so redo() should go to index 2
   int modification_index = 1;
   
   // Call redo_gemmi()
   const gemmi::Structure* result = coot::molecule_t::modification_info_t::redo_gemmi(save_info, modification_index);
   
   // Verify that redo_gemmi() returned a non-null structure
   EXPECT_NE(result, nullptr);
   // Verify that modification_index was incremented to 2
   EXPECT_EQ(modification_index, 2);
   
   // Test edge case: redo when already at end of history
   modification_index = save_info.size() - 1;  // Set to last index (2)
   const gemmi::Structure* result2 = coot::molecule_t::modification_info_t::redo_gemmi(save_info, modification_index);
   
   // redo_gemmi() should return null when already at end of history
   EXPECT_EQ(result2, nullptr);
   // modification_index should remain unchanged at end of history
   EXPECT_EQ(modification_index, save_info.size() - 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}