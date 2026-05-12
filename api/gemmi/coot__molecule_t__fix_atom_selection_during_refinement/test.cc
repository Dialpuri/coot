#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(FixAtomSelectionTest, Basic) {
   // Load test structure
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   coot::molecule_t mol(st);
   
   // Test case 1: CID //A/10
   {
      std::string cid1 = "//A/10";
      std::vector<coot::atom_spec_t> fixed_before = mol.get_fixed_atoms();
      EXPECT_EQ(fixed_before.size(), 0);
      
      mol.fix_atom_selection_during_refinement_gemmi(cid1);
      
      std::vector<coot::atom_spec_t> fixed_after = mol.get_fixed_atoms();
      EXPECT_EQ(fixed_after.size(), 0);
   }
   
   // Test case 2: CID //A/15
   {
      std::string cid2 = "//A/15";
      mol.fix_atom_selection_during_refinement_gemmi(cid2);
      
      std::vector<coot::atom_spec_t> fixed_after2 = mol.get_fixed_atoms();
      EXPECT_EQ(fixed_after2.size(), 0);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}