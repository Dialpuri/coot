#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace coot {

TEST(AddNTerminalResidueBondsAndAnglesToHydrogens, MET_N_terminal) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
  gemmi::Model &model = st.models[0];

  // Find MET A/1
  gemmi::Residue *met_res = nullptr;
  for (gemmi::Chain &chain : model.chains) {
    if (chain.name == "A") {
      for (gemmi::Residue &res : chain.residues) {
        if (res.name == "MET" && res.seqid.num.value == 1) {
          met_res = &res;
          break;
        }
      }
    }
  }
  ASSERT_NE(met_res, nullptr);

  coot::restraints_container_t::restraint_counts_t rc_counts =
      add_N_terminal_residue_bonds_and_angles_to_hydrogens_gemmi(met_res);

  EXPECT_EQ(rc_counts.n_bond_restraints, 2);
  EXPECT_EQ(rc_counts.n_angle_restraints, 2);
}

TEST(AddNTerminalResidueBondsAndAnglesToHydrogens, VAL_no_hydrogens) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
  gemmi::Model &model = st.models[0];

  // Find VAL A/7
  gemmi::Residue *val_res = nullptr;
  for (gemmi::Chain &chain : model.chains) {
    if (chain.name == "A") {
      for (gemmi::Residue &res : chain.residues) {
        if (res.name == "VAL" && res.seqid.num.value == 7) {
          val_res = &res;
          break;
        }
      }
    }
  }
  ASSERT_NE(val_res, nullptr);

  coot::restraints_container_t::restraint_counts_t rc_counts =
      add_N_terminal_residue_bonds_and_angles_to_hydrogens_gemmi(val_res);

  EXPECT_EQ(rc_counts.n_bond_restraints, 0);
  EXPECT_EQ(rc_counts.n_angle_restraints, 0);
}

TEST(AddNTerminalResidueBondsAndAnglesToHydrogens, LYS_no_hydrogens) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
  gemmi::Model &model = st.models[0];

  // Find LYS A/9
  gemmi::Residue *lys_res = nullptr;
  for (gemmi::Chain &chain : model.chains) {
    if (chain.name == "A") {
      for (gemmi::Residue &res : chain.residues) {
        if (res.name == "LYS" && res.seqid.num.value == 9) {
          lys_res = &res;
          break;
        }
      }
    }
  }
  ASSERT_NE(lys_res, nullptr);

  coot::restraints_container_t::restraint_counts_t rc_counts =
      add_N_terminal_residue_bonds_and_angles_to_hydrogens_gemmi(lys_res);

  EXPECT_EQ(rc_counts.n_bond_restraints, 0);
  EXPECT_EQ(rc_counts.n_angle_restraints, 0);
}

} // namespace coot

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
