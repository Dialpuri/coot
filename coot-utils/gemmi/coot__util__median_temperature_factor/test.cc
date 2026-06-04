#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, median_temperature_factor) {
  gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  // Collect all atoms from the structure
  std::vector<gemmi::Atom const*> atom_selection;
  for (gemmi::Model& model : st.models) {
    for (gemmi::Chain& chain : model.chains) {
      for (gemmi::Residue& residue : chain.residues) {
        for (gemmi::Atom& atom : residue.atoms) {
          atom_selection.push_back(&atom);
        }
      }
    }
  }

  int n_atoms = static_cast<int>(atom_selection.size());
  // The original test expects 2107 atoms
  EXPECT_EQ(n_atoms, 2107);

  // Case 1: no cutoffs applied (default behaviour)
  {
    float low_cutoff = 2.0f;
    float high_cutoff = 2222.2f;
    bool apply_low = false;
    bool apply_high = false;
    float median = coot::util::median_temperature_factor_gemmi(
        atom_selection, low_cutoff, high_cutoff, apply_low, apply_high);
    EXPECT_NEAR(median, 28.19, 1e-4);
  }

  // Case 2: apply low cutoff (exclude atoms with b <= 2.0)
  {
    float low_cutoff = 2.0f;
    float high_cutoff = 100.0f;
    bool apply_low = true;
    bool apply_high = false;
    float median = coot::util::median_temperature_factor_gemmi(
        atom_selection, low_cutoff, high_cutoff, apply_low, apply_high);
    EXPECT_NEAR(median, 28.19, 1e-4);
  }

  // Case 3: apply both cutoffs (exclude atoms with b <= 0 or b >= 100)
  {
    float low_cutoff = 0.0f;
    float high_cutoff = 100.0f;
    bool apply_low = true;
    bool apply_high = true;
    float median = coot::util::median_temperature_factor_gemmi(
        atom_selection, low_cutoff, high_cutoff, apply_low, apply_high);
    EXPECT_NEAR(median, 103.14, 1e-4);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}