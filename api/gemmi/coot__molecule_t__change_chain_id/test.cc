#define _FNO_ACCESS_CONTROL
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ChangeChainID) {
   // Load the test PDB file
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Test case 1: Simple chain ID change (A -> X) - should succeed
   {
      auto result = coot::molecule_t::change_chain_id_gemmi(st, "A", "X", false, 0, 0);
      EXPECT_EQ(result.first, 1);
      // Status 1 means success
   }

   // Reload the structure (since we changed it)
   st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Test case 2: Chain ID change with residue range (B -> Y, residues 1-10)
   {
      auto result = coot::molecule_t::change_chain_id_gemmi(st, "B", "Y", true, 1, 10);
      EXPECT_EQ(result.first, 1);
   }

   // Reload the structure
   st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Test case 3: Conflict case - target chain already exists (A -> B)
   {
      auto result = coot::molecule_t::change_chain_id_gemmi(st, "A", "B", false, 0, 0);
      EXPECT_EQ(result.first, 0);
      // Status 0 means failure due to conflict
      EXPECT_NE(result.second.find("CONFLICT"), std::string::npos);
   }

   // Test case 4: Non-existent source chain (Z -> W)
   {
      auto result = coot::molecule_t::change_chain_id_gemmi(st, "Z", "W", false, 0, 0);
      EXPECT_EQ(result.first, 0);
   }

   // Test case 5: Same chain ID (A -> A) - should do nothing but succeed
   {
      auto result = coot::molecule_t::change_chain_id_gemmi(st, "A", "A", false, 0, 0);
      EXPECT_EQ(result.first, 0);
   }
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}