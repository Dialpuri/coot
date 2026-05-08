#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetBondsMesh) {
   coot::molecule_t mol("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
   
   struct test_case_t {
      std::string mode;
      bool against_a_dark_background;
      float bonds_width;
      float atom_radius_to_bond_width_ratio;
      int smoothness_factor;
      bool draw_hydrogen_atoms_flag;
      size_t expected_vertices;
      size_t expected_triangles;
      size_t expected_colour_map_size;
   };
   
   std::vector<test_case_t> test_cases = {
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.0f, 2, true, 720066, 1044128, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.0f, 1, true, 303226, 385104, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.0f, 3, true, 1977466, 3183936, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.0f, 4, true, 6122946, 10750592, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.0f, 2, false, 381078, 557376, 52},
      {"CA+LIGANDS", true, 0.3f, 1.0f, 2, true, 48708, 57024, 0},
      {"VDW-BALLS", true, 0.3f, 1.0f, 2, true, 316720, 547840, 0},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.2f, 1.0f, 2, true, 591622, 770144, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.4f, 1.0f, 2, true, 591622, 770144, 52},
      {"COLOUR-BY-CHAIN-AND-DICTIONARY", true, 0.3f, 1.5f, 2, true, 720066, 1044128, 52},
   };
   
   for (size_t i = 0; i < test_cases.size(); ++i) {
      const auto& tc = test_cases[i];
      
      coot::simple_mesh_t mesh = mol.get_bonds_mesh_gemmi(tc.mode, tc.against_a_dark_background, 
                                                         tc.bonds_width, tc.atom_radius_to_bond_width_ratio,
                                                         tc.smoothness_factor, tc.draw_hydrogen_atoms_flag, false);
      
      EXPECT_EQ(mesh.vertices.size(), tc.expected_vertices);
      EXPECT_EQ(mesh.triangles.size(), tc.expected_triangles);
      EXPECT_EQ(mesh.status, 1);
      EXPECT_EQ(mesh.colour_index_to_colour_map.size(), tc.expected_colour_map_size);
   }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}