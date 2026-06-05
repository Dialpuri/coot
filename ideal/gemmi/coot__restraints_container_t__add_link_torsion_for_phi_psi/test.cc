#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"
#include "clipper/core/xmap.h"

TEST(OracleTest, add_link_torsion_for_phi_psi) {
  gemmi::Structure st = gemmi::read_pdb_file(
      "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

  gemmi::Model &model = st.models[0];

  // Build cras from the model
  std::vector<gemmi::CRA> cras;
  for (gemmi::Chain &chain : model.chains) {
    for (gemmi::Residue &residue : chain.residues) {
      for (gemmi::Atom &atom : residue.atoms) {
        cras.push_back({&chain, &residue, &atom});
      }
    }
  }

  coot::protein_geometry pg;
  pg.init_standard();

  clipper::Xmap<float> dummy_xmap;

  // Find residues 5 and 6 in chain A
  gemmi::CRA cra_res1 = {nullptr, nullptr, nullptr};
  gemmi::CRA cra_res2 = {nullptr, nullptr, nullptr};

  for (gemmi::Chain &chain : model.chains) {
    if (chain.name == "A") {
      for (gemmi::Residue &residue : chain.residues) {
        if (residue.seqid.num.value == 5) {
          cra_res1 = {&chain, &residue, nullptr};
        }
        if (residue.seqid.num.value == 6) {
          cra_res2 = {&chain, &residue, nullptr};
        }
      }
    }
  }

  EXPECT_NE(cra_res1.residue, nullptr);
  EXPECT_NE(cra_res2.residue, nullptr);

  std::string link_type = "TRANS";
  int phi_psi_type = coot::restraints_container_t::LINK_TORSION_ALPHA_HELIX;
  short int is_fixed_first = 0;
  short int is_fixed_second = 0;

  // Case 1: TRANS + ALPHA_HELIX
  {
    std::vector<coot::simple_restraint> restraints_vec;
    int n_added = coot::add_link_torsion_for_phi_psi_gemmi(
        restraints_vec, cras,
        link_type, phi_psi_type, cra_res1, cra_res2,
        is_fixed_first, is_fixed_second, pg);

    EXPECT_EQ(n_added, 2);
  }

  // Case 2: BETA_STRAND
  {
    int phi_psi_type2 = coot::restraints_container_t::LINK_TORSION_BETA_STRAND;

    std::vector<coot::simple_restraint> restraints_vec;
    int n_added2 = coot::add_link_torsion_for_phi_psi_gemmi(
        restraints_vec, cras,
        link_type, phi_psi_type2, cra_res1, cra_res2,
        is_fixed_first, is_fixed_second, pg);

    EXPECT_EQ(n_added2, 2);
  }

  // Case 3: invalid link_type — should add 0
  {
    std::string link_type3 = "INVALID_TYPE";

    std::vector<coot::simple_restraint> restraints_vec;
    int n_added3 = coot::add_link_torsion_for_phi_psi_gemmi(
        restraints_vec, cras,
        link_type3, phi_psi_type, cra_res1, cra_res2,
        is_fixed_first, is_fixed_second, pg);

    EXPECT_EQ(n_added3, 0);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}