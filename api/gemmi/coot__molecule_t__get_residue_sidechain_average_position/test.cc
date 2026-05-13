#define _NO_GTKMM
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_residue_sidechain_average_position) {
   // Setup - load the structure directly since we don't have molecules_container_t
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   gemmi::Model& model = st.models[0];

   // Test 1: Valid residue with sidechain (MET 1, chain A)
   {  // Scope for test 1
      std::string cid1 = "//A/1";
      std::vector<double> result1 = coot::molecule_t::get_residue_sidechain_average_position_gemmi(cid1, model);
      EXPECT_EQ(result1.size(), 3u);
      EXPECT_NEAR(result1[0], 13.9735, 0.0001);
      EXPECT_NEAR(result1[1], 6.531, 0.0001);
      EXPECT_NEAR(result1[2], 74.5918, 0.0001);
   }

   // Test 2: Another valid residue (ALA 31, chain A)
   {  // Scope for test 2
      std::string cid2 = "//A/31";
      std::vector<double> result2 = coot::molecule_t::get_residue_sidechain_average_position_gemmi(cid2, model);
      EXPECT_EQ(result2.size(), 3u);
      EXPECT_NEAR(result2[0], 25.158, 0.0001);
      EXPECT_NEAR(result2[1], 4.628, 0.0001);
      EXPECT_NEAR(result2[2], 67.262, 0.0001);
   }

   // Test 3: Glycine (no sidechain)
   {  // Scope for test 3
      std::string cid3 = "//A/11";
      std::vector<double> result3 = coot::molecule_t::get_residue_sidechain_average_position_gemmi(cid3, model);
      EXPECT_EQ(result3.size(), 0u);
   }

   // Test 4: Invalid CID (nonexistent residue)
   {  // Scope for test 4
      std::string cid4 = "//X/9999";
      std::vector<double> result4 = coot::molecule_t::get_residue_sidechain_average_position_gemmi(cid4, model);
      EXPECT_EQ(result4.size(), 0u);
   }

   // Test 5: Residue with large sidechain (TRP 167, chain A)
   {  // Scope for test 5
      std::string cid5 = "//A/167";
      std::vector<double> result5 = coot::molecule_t::get_residue_sidechain_average_position_gemmi(cid5, model);
      EXPECT_EQ(result5.size(), 3u);
      EXPECT_NEAR(result5[0], 30.599, 0.0001);
      EXPECT_NEAR(result5[1], 15.126, 0.0001);
      EXPECT_NEAR(result5[2], 46.75, 0.0001);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}