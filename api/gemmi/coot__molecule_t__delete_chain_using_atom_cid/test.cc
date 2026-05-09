#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper function to count chains in a structure
static int count_chains(const gemmi::Structure& st) {
   int nChains = 0;
   for (const auto& model : st.models) {
      nChains += (int)model.chains.size();
   }
   return nChains;
}

TEST(OracleTest, delete_chain_using_atom_cid) {
   // Load initial PDB (2 chains: A and B)
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Verify initial state: 2 chains
   int nChainsBefore = count_chains(st);
   EXPECT_EQ(nChainsBefore, 2);
   
   // Test 1: Delete chain 'A' using atom CID in chain A
   {
      std::string cid1 = "//A/1/CA";
      int result1 = coot::delete_chain_using_atom_cid_gemmi(st, cid1);
      EXPECT_EQ(result1, 1);  // operation succeeded
      
      int nChainsAfter1 = count_chains(st);
      EXPECT_EQ(nChainsAfter1, 1);  // 1 chain left (B)
   }
   
   // Load fresh structure for second test
   gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Test 2: Delete chain 'B' (water chain) using atom CID in chain B
   {
      std::string cid2 = "//B/1/O";
      int result2 = coot::delete_chain_using_atom_cid_gemmi(st2, cid2);
      EXPECT_EQ(result2, 1);  // operation succeeded
      
      int nChainsAfter2 = count_chains(st2);
      EXPECT_EQ(nChainsAfter2, 1);  // 1 chain left (A)
   }
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}