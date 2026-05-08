#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, cid_to_residue) {
   // Setup: load molecule
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   gemmi::Model& model = st.models[0];

   // Case 1: Valid CID for chain A, residue 10 (ILE)
   {
      std::string cid_valid = "//A/10";
      gemmi::Residue* res = coot::cid_to_residue_gemmi(cid_valid, model);
      EXPECT_TRUE(res != nullptr);
      if (res) {
         EXPECT_EQ(res->seqid.num.value, 10);
         EXPECT_STRNE(res->name.c_str(), "");
      }
   }

   // Case 2: Invalid CID (non-existent residue)
   {
      std::string cid_invalid = "//A/999";
      gemmi::Residue* res = coot::cid_to_residue_gemmi(cid_invalid, model);
      EXPECT_TRUE(res == nullptr);
   }

   // Case 3: Chain B, residue 1 (HOH)
   {
      std::string cid_b = "//B/1";
      gemmi::Residue* res = coot::cid_to_residue_gemmi(cid_b, model);
      EXPECT_TRUE(res != nullptr);
      if (res) {
         EXPECT_EQ(res->seqid.num.value, 1);
         EXPECT_STRNE(res->name.c_str(), "");
      }
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}