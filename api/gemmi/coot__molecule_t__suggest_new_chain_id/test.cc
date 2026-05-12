#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, suggest_new_chain_id) {
   // Test Case 1: Chain 'A' exists, should suggest next available
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      std::string input1 = "A";
      std::string output1 = coot::molecule_t_suggest_new_chain_id_gemmi(st, input1);
      
      // 'A' exists, so it should suggest next available (not 'A')
      EXPECT_FALSE(output1.empty());
      EXPECT_NE(output1, "A");
   }
   
   // Test Case 2: Chain 'B' exists (water), should suggest next available
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      std::string input2 = "B";
      std::string output2 = coot::molecule_t_suggest_new_chain_id_gemmi(st, input2);
      
      // 'B' exists (water), so should suggest next available
      EXPECT_FALSE(output2.empty());
      EXPECT_NE(output2, "B");
   }
   
   // Test Case 3: Chain 'X' doesn't exist, should return 'X' directly
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      std::string input3 = "X";
      std::string output3 = coot::molecule_t_suggest_new_chain_id_gemmi(st, input3);
      
      // 'X' doesn't exist, should return 'X' directly
      EXPECT_EQ(output3, "X");
   }
   
   // Test Case 4: Multi-char chain ID (test postfix logic)
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      std::string input4 = "AB";
      std::string output4 = coot::molecule_t_suggest_new_chain_id_gemmi(st, input4);
      
      // 'AB' doesn't exist, should return 'AB' directly
      EXPECT_EQ(output4, "AB");
   }
   
   // Test Case 5: Chain 'C' doesn't exist, should return 'C' directly
   {
      gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
      
      std::string input5 = "C";
      std::string output5 = coot::molecule_t_suggest_new_chain_id_gemmi(st, input5);
      
      // 'C' doesn't exist, should return 'C' directly
      EXPECT_EQ(output5, "C");
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}