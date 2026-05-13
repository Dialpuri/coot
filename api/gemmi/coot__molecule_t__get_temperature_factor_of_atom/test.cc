#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include <gemmi/pdb.hpp>
#include "function.hh"

TEST(OracleTest, GetTemperatureFactorOfAtom) {
   // Load the test structure
   std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
   gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

   // Test 1: Valid atom CID //A/10/CA -> expected b_factor = 41.55
   {
      std::string atom_cid1 = "//A/10/CA";
      float b1 = coot::molecule_t_get_temperature_factor_of_atom_gemmi(atom_cid1, st);
      EXPECT_FLOAT_EQ(b1, 41.55f);
   }

   // Test 2: Valid atom CID //A/1/N -> expected b_factor = 65.05
   {
      std::string atom_cid2 = "//A/1/N";
      float b2 = coot::molecule_t_get_temperature_factor_of_atom_gemmi(atom_cid2, st);
      EXPECT_FLOAT_EQ(b2, 65.05f);
   }

   // Test 3: Invalid atom CID //X/9999/ZZZ -> expected b_factor = -1.1
   {
      std::string atom_cid_invalid = "//X/9999/ZZZ";
      float b_invalid = coot::molecule_t_get_temperature_factor_of_atom_gemmi(atom_cid_invalid, st);
      EXPECT_FLOAT_EQ(b_invalid, -1.1f);
   }
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}