#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residue_is_nucleic_acid) {
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Test case 1: Protein residue (should be false)
   {
      std::string cid1 = "//A/10";  // ILE 10
      bool result1 = coot::molecule_t::residue_is_nucleic_acid_gemmi(cid1, st);
      EXPECT_FALSE(result1);
   }
   
   // Test case 2: Another protein residue (should be false)
   {
      std::string cid2 = "//A/20";  // LYS 20
      bool result2 = coot::molecule_t::residue_is_nucleic_acid_gemmi(cid2, st);
      EXPECT_FALSE(result2);
   }
   
   // Test case 3: Non-existent residue (should be false)
   {
      std::string cid3 = "//A/9999";  // Non-existent
      bool result3 = coot::molecule_t::residue_is_nucleic_acid_gemmi(cid3, st);
      EXPECT_FALSE(result3);
   }
   
   // Test case 4: Different chain (should be false)
   {
      std::string cid4 = "//B/10";  // Chain B (water)
      bool result4 = coot::molecule_t::residue_is_nucleic_acid_gemmi(cid4, st);
      EXPECT_FALSE(result4);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}