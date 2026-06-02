#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string strip_atom_name(const std::string& name) {
    std::string trimmed = name;
    trimmed.erase(0, trimmed.find_first_not_of(' '));
    if (!trimmed.empty())
        trimmed.erase(trimmed.find_last_not_of(' ') + 1);
    return trimmed;
}

TEST(OracleTest, new_positions_for_residue_atoms) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model& model = st.models[0];

  // Case 1: Move one atom (CA) in residue //A/1 (MET)
  {
    std::string residue_cid = "//A/1";

    // Get original CA position
    gemmi::Residue* res = coot::cid_to_residue_gemmi(residue_cid, model);

    double orig_x = 0, orig_y = 0, orig_z = 0;
    std::string ca_name;
    if (res) {
      for (auto& atom : res->atoms) {
        std::string aname = atom.name;
        if (strip_atom_name(aname) == "CA") {
          orig_x = atom.pos.x;
          orig_y = atom.pos.y;
          orig_z = atom.pos.z;
          ca_name = aname;
          break;
        }
      }
    }

    EXPECT_NEAR(orig_x, 13.401, 1e-4);
    EXPECT_NEAR(orig_y, 3.78, 1e-4);
    EXPECT_NEAR(orig_z, 74.703, 1e-4);

    std::vector<coot::api::moved_atom_t> moved_atoms;
    coot::api::moved_atom_t ca_move(ca_name, "", orig_x + 1.0f, orig_y + 1.0f, orig_z + 1.0f);
    moved_atoms.push_back(ca_move);

    EXPECT_EQ(ca_move.atom_name, "CA");
    EXPECT_NEAR(ca_move.x, 14.401, 1e-4);
    EXPECT_NEAR(ca_move.y, 4.78, 1e-4);
    EXPECT_NEAR(ca_move.z, 75.703, 1e-4);

    int result = coot::new_positions_for_residue_atoms_gemmi(residue_cid, moved_atoms, model);
    EXPECT_EQ(result, 1);

    // Verify new position
    if (res) {
      for (auto& atom : res->atoms) {
        std::string aname = atom.name;
        if (strip_atom_name(aname) == "CA") {
          EXPECT_NEAR(atom.pos.x, 14.401, 1e-4);
          EXPECT_NEAR(atom.pos.y, 4.78, 1e-4);
          EXPECT_NEAR(atom.pos.z, 75.703, 1e-4);
          break;
        }
      }
    }
  }

  // Case 2: Invalid residue CID
  {
    std::string residue_cid = "//A/9999";

    std::vector<coot::api::moved_atom_t> moved_atoms;
    coot::api::moved_atom_t ca_move("CA", "", 0.0f, 0.0f, 0.0f);
    moved_atoms.push_back(ca_move);

    int result = coot::new_positions_for_residue_atoms_gemmi(residue_cid, moved_atoms, model);
    EXPECT_EQ(result, -1);
  }

  // Case 3: Empty moved_atoms vector (residue exists but nothing to move)
  {
    std::string residue_cid = "//A/1";
    std::vector<coot::api::moved_atom_t> moved_atoms; // empty

    int result = coot::new_positions_for_residue_atoms_gemmi(residue_cid, moved_atoms, model);
    EXPECT_EQ(result, 0);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}