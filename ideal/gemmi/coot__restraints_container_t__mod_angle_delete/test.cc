#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static int count_angle_restraints_for_residue(
    const restraints_container_t_gemmi &rc, gemmi::Residue *residue_p) {
  int count = 0;
  for (const auto &r : rc.restraints_vec) {
    if (r.size() >= 4 && r[0] == ANGLE_RESTRAINT_TYPE) {
      size_t i1 = static_cast<size_t>(r[1]);
      size_t i2 = static_cast<size_t>(r[2]);
      if (i1 < rc.atom_residues.size() && i2 < rc.atom_residues.size() &&
          rc.atom_residues[i1] == residue_p &&
          rc.atom_residues[i2] == residue_p) {
        count++;
      }
    }
  }
  return count;
}

static std::string trim_atom_name(std::string name) {
  auto lo = name.find_first_not_of(' ');
  if (lo == std::string::npos) return "";
  auto hi = name.find_last_not_of(' ');
  return name.substr(lo, hi - lo + 1);
}

TEST(OracleTest, mod_angle_delete) {
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
  gemmi::Model &model = st.models[0];

  // Find chain A, residue 10
  gemmi::Residue *target_residue = nullptr;
  for (auto &ch : model.chains) {
    for (auto &res : ch.residues) {
      if (ch.name == "A" && res.seqid.num.value == 10) {
        target_residue = &res;
      }
    }
  }
  ASSERT_NE(target_residue, nullptr);
  ASSERT_GT(target_residue->atoms.size(), 4u);

  // Build restraints_container
  restraints_container_t_gemmi rc;
  for (size_t i = 0; i < target_residue->atoms.size(); ++i) {
    rc.atom_residues.push_back(target_residue);
    rc.atom_names.push_back(trim_atom_name(target_residue->atoms[i].name));
  }

  // Atom indices for this residue (N=0, CA=1, C=2, O=3, CB=4, CG1=5, CG2=6, CD1=7)
  int n_idx = 0, ca_idx = 1, c_idx = 2, o_idx = 3, cb_idx = 4;

  // Add 10 angle restraints that do NOT have (N, CA) as atom_1/atom_2
  for (int k = 0; k < 10; ++k) {
    std::vector<int> ang;
    ang.push_back(ANGLE_RESTRAINT_TYPE);
    ang.push_back(o_idx); // O as atom_1 (not "N")
    ang.push_back(cb_idx); // CB as atom_2 (not "CA")
    ang.push_back(c_idx);
    rc.restraints_vec.push_back(ang);
  }

  // === Frozen oracle test cases ===

  // Case 1: mod_angle("DELETE", "N", "CA", "CA") — no matching restraint → stays 10
  {
    int angle_count_before = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_before, 10);

    chem_mod_angle_gemmi mod_angle("DELETE", "N", "CA", "CA", 0, 0);
    rc.mod_angle_delete_gemmi(mod_angle, target_residue);

    int angle_count_after = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_after, 10);
  }

  // Case 2: Non-matching atom_id_1 ("OX1" → no match)
  {
    chem_mod_angle_gemmi mod_angle2("DELETE", "OX1", "CA", "CA", 0, 0);
    int angle_count_before2 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_before2, 10);

    rc.mod_angle_delete_gemmi(mod_angle2, target_residue);

    int angle_count_after2 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_after2, 10);
  }

  // === Additional: verify deletion actually works ===

  // Case 3: Add N-CA-C restraint, then delete with atom_id_3=="CA" (buggy match)
  {
    std::vector<int> ang;
    ang.push_back(ANGLE_RESTRAINT_TYPE);
    ang.push_back(n_idx);
    ang.push_back(ca_idx);
    ang.push_back(c_idx);
    rc.restraints_vec.push_back(ang);

    int angle_count_before3 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_before3, 11);

    // Buggy: name_2("CA") == atom_id_3("CA") → TRUE → deletes
    chem_mod_angle_gemmi mod_angle3("DELETE", "N", "CA", "CA", 0, 0);
    rc.mod_angle_delete_gemmi(mod_angle3, target_residue);

    int angle_count_after3 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_after3, 10);
  }

  // Case 4: atom_id_3 != atom_id_2 → buggy check fails → no delete
  {
    std::vector<int> ang;
    ang.push_back(ANGLE_RESTRAINT_TYPE);
    ang.push_back(n_idx);
    ang.push_back(ca_idx);
    ang.push_back(c_idx);
    rc.restraints_vec.push_back(ang);

    // Buggy: name_2("CA") == atom_id_3("C") → FALSE → no delete
    chem_mod_angle_gemmi mod_angle4("DELETE", "N", "CA", "C", 0, 0);
    int angle_count_before4 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_before4, 11);

    rc.mod_angle_delete_gemmi(mod_angle4, target_residue);

    int angle_count_after4 = count_angle_restraints_for_residue(rc, target_residue);
    EXPECT_EQ(angle_count_after4, 11);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}