#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cid_to_residue_spec) {
   using namespace coot;
   
   // Load structure using gemmi
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   {  // Test with a valid CID that should match a residue in chain A
      std::string cid_valid = "//A/10";
      std::pair<bool, residue_spec_t> result_valid = cid_to_residue_spec_gemmi(st, cid_valid);
      EXPECT_TRUE(result_valid.first);
      EXPECT_EQ(result_valid.second.chain_id, "A");
      EXPECT_EQ(result_valid.second.res_no, 10);
      EXPECT_EQ(result_valid.second.ins_code, "");
   }

   {  // Test with an invalid CID
      std::string cid_invalid = "//X/999";
      std::pair<bool, residue_spec_t> result_invalid = cid_to_residue_spec_gemmi(st, cid_invalid);
      EXPECT_FALSE(result_invalid.first);
   }
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}