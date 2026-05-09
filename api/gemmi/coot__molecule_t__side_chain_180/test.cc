#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, side_chain_180) {
   // Load the structure using gemmi
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Test 1: ALA 1 (small side chain)
   {
      int status = coot::molecule_t::side_chain_180_gemmi(st, "A", 1, "");
      EXPECT_EQ(status, 1);
   }

   // Test 2: GLU 2 (larger side chain with multiple chi angles)
   {
      int status = coot::molecule_t::side_chain_180_gemmi(st, "A", 2, "");
      EXPECT_EQ(status, 1);
   }

   // Test 3: ARG 22 (large side chain with 5 chi angles)
   {
      int status = coot::molecule_t::side_chain_180_gemmi(st, "A", 22, "");
      EXPECT_EQ(status, 1);
   }

   // Test 4: GLY 5 (no side chain)
   {
      int status = coot::molecule_t::side_chain_180_gemmi(st, "A", 5, "");
      EXPECT_EQ(status, 1);
   }

   // Test 5: SER 10 (hydroxyl side chain)
   {
      int status = coot::molecule_t::side_chain_180_gemmi(st, "A", 10, "");
      EXPECT_EQ(status, 1);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}