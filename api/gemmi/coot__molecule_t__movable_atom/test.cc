#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MovableAtom) {
   // Load a PDB file with explicit hydrogens to have zero-occupancy atoms
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
   
   ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file or no models found";
   
   // Get the first atom (should be non-zero occupancy)
   const gemmi::Atom *atom = nullptr;
   for (auto& model : st.models)
      for (auto& chain : model.chains)
         for (auto& res : chain.residues)
            for (auto& a : res.atoms) {
               atom = &a;
               goto done;
            }
   done:
   
   ASSERT_NE(atom, nullptr) << "No atoms found in structure";
   
   // Test 1: non-zero occupancy atom with replace_coords_with_zero_occ_flag = false
   bool result1 = coot::movable_atom_gemmi(atom, false);
   EXPECT_TRUE(result1);
   
   // Test 2: non-zero occupancy atom with replace_coords_with_zero_occ_flag = true
   bool result2 = coot::movable_atom_gemmi(atom, true);
   EXPECT_TRUE(result2);
   
   // Test 3: zero occupancy atom with replace_coords_with_zero_occ_flag = false
   const gemmi::Atom *zero_occ_atom = nullptr;
   for (auto& model : st.models)
      for (auto& chain : model.chains)
         for (auto& res : chain.residues)
            for (auto& a : res.atoms)
               if (a.occ < 0.0001 && a.occ > -0.0001) {
                  zero_occ_atom = &a;
                  goto done2;
               }
   done2:
   
   if (zero_occ_atom) {
      bool result3 = coot::movable_atom_gemmi(zero_occ_atom, false);
      EXPECT_FALSE(result3);
      
      bool result4 = coot::movable_atom_gemmi(zero_occ_atom, true);
      EXPECT_TRUE(result4);
   }
   
   // Test 5: null atom (error case)
   bool result5 = coot::movable_atom_gemmi(nullptr, false);
   EXPECT_FALSE(result5);
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}