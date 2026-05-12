#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include <tuple>
#include <cmath>
#include "function.hh"

TEST(OracleTest, SetResidueToRotamerMoveAtoms) {
   using namespace coot;
   
   // Load molecule via gemmi
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
   
   // Get residue 2 from chain A (GLU in example.pdb)
   gemmi::Model& model = st.models[0];
   gemmi::Chain& chain = model.chains[0];
   gemmi::Residue& res = chain.residues[1];  // index 1 = residue 2
   
   ASSERT_EQ(res.name, "GLU");
   
   // Store original coordinates for comparison
   std::map<std::string, std::tuple<float, float, float>> original_coords;
   for (const gemmi::Atom& atom : res.atoms) {
      original_coords[atom.name] = std::make_tuple(atom.pos.x, atom.pos.y, atom.pos.z);
   }
   
   // Create a moving residue by copying res, which ensures atom names match
   gemmi::Residue moving_res = res;
   // Modify some atom coordinates in moving_res
   for (gemmi::Atom& atom : moving_res.atoms) {
      atom.pos.x += 1.0;
      atom.pos.y += 1.0;
      atom.pos.z += 1.0;
   }
   
   // Record before state
   float before_x = res.atoms[0].pos.x;
   float before_y = res.atoms[0].pos.y;
   float before_z = res.atoms[0].pos.z;
   
   // Call the function
   int result = molecule_t::set_residue_to_rotamer_move_atoms_gemmi(res, moving_res);
   
   // Verify return value indicates success
   EXPECT_EQ(result, 1);
   
   // Record after state
   int n_ref_atoms = res.atoms.size();
   int n_after = res.atoms.size();
   EXPECT_EQ(n_after, n_ref_atoms);
   
   // Verify that coordinates have been updated
   // All atoms should have changed since we moved all atoms by 1.0 in each direction
   int atoms_changed = 0;
   for (size_t i = 0; i < res.atoms.size(); ++i) {
      std::string atom_name = res.atoms[i].name;
      auto it = original_coords.find(atom_name);
      if (it != original_coords.end()) {
         float orig_x = std::get<0>(it->second);
         float orig_y = std::get<1>(it->second);
         float orig_z = std::get<2>(it->second);
         float dx = std::abs(res.atoms[i].pos.x - orig_x - 1.0f);
         float dy = std::abs(res.atoms[i].pos.y - orig_y - 1.0f);
         float dz = std::abs(res.atoms[i].pos.z - orig_z - 1.0f);
         if (dx < 0.01f && dy < 0.01f && dz < 0.01f) {
            atoms_changed++;
         }
      }
   }
   
   // All atoms should have been moved (copied from moving_res)
   EXPECT_EQ(atoms_changed, (int)res.atoms.size());
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}