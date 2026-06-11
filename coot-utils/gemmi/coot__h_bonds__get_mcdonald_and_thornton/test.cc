#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

// Helper: build a vector of CRA from selection criteria on a gemmi Model
static std::vector<gemmi::CRA> select_atoms(
  gemmi::Model& model,
  const std::string& chain_id,
  int resno_min,
  const std::string& icode_min,
  int resno_max,
  const std::string& icode_max,
  bool select_all = false) {

  std::vector<gemmi::CRA> result;

  if (select_all) {
    for (auto& chain : model.chains) {
      for (auto& residue : chain.residues) {
        for (auto& atom : residue.atoms) {
          gemmi::CRA cra{&chain, &residue, &atom};
          result.push_back(cra);
        }
      }
    }
    return result;
  }

  // Normalize icode: empty string in MMDB means ' ' in gemmi
  auto norm = [](const std::string& ic) -> std::string {
    return ic.empty() ? std::string(" ") : ic;
  };
  std::string norm_icode_min = norm(icode_min);
  std::string norm_icode_max = norm(icode_max);

  for (auto& chain : model.chains) {
    if (chain_id != "*" && chain.name != chain_id) continue;
    for (auto& residue : chain.residues) {
      int seq = residue.seqid.num.value;
      std::string icode(1, residue.seqid.icode);

      if (seq >= resno_min && seq <= resno_max &&
          icode == norm_icode_min && icode == norm_icode_max) {
        for (auto& atom : residue.atoms) {
          gemmi::CRA cra{&chain, &residue, &atom};
          result.push_back(cra);
        }
      }
    }
  }
  return result;
}

// Helper: build atom* -> chain_name, resno, atom_name lookup (gemmi atoms have no parent pointer)
static void build_atom_info(
    gemmi::Model& model,
    std::map<const gemmi::Atom*, std::string>& atom_chain,
    std::map<const gemmi::Atom*, int>& atom_resno) {
  for (auto& chain : model.chains) {
    for (auto& residue : chain.residues) {
      for (auto& atom : residue.atoms) {
        atom_chain[&atom] = chain.name;
        atom_resno[&atom] = residue.seqid.num.value;
      }
    }
  }
}

TEST(OracleTest, get_mcdonald_and_thornton) {
  gemmi::Structure st = gemmi::read_pdb_file(
    "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

  gemmi::Model& model = st.models[0];

  coot::protein_geometry geom;
  geom.init_standard();

  // Build atom info lookup (gemmi atoms have no parent pointer)
  std::map<const gemmi::Atom*, std::string> atom_chain;
  std::map<const gemmi::Atom*, int> atom_resno;
  build_atom_info(model, atom_chain, atom_resno);

  // Case 1: Use a polar residue (SER 46) as "ligand" selection, all atoms as environment
  {
    auto sel_all = select_atoms(model, "*", 0, "", 99999, "", true);
    auto sel_lig = select_atoms(model, "A", 46, "", 46, "");

    double max_dist = 3.8;

    auto hbonds = coot::get_mcdonald_and_thornton_gemmi(
        sel_lig, sel_all, model, geom, 0, max_dist);

    EXPECT_EQ(hbonds.size(), 2u);

    if (hbonds.size() >= 1) {
      EXPECT_NEAR(hbonds[0].dist, 2.50604, 1e-4);
      EXPECT_NEAR(hbonds[0].angle_1, 133.445, 0.001);
      EXPECT_NEAR(hbonds[0].angle_2, -1.0, 1e-4);
      EXPECT_NEAR(hbonds[0].angle_3, 141.798, 0.001);
      EXPECT_EQ(hbonds[0].ligand_atom_is_donor, false);
      EXPECT_EQ(hbonds[0].bond_has_hydrogen_flag, true);

      // donor chain A, res 50, atom name contains "N"
      if (hbonds[0].donor) {
        EXPECT_EQ(atom_chain[hbonds[0].donor], "A");
        EXPECT_EQ(atom_resno[hbonds[0].donor], 50);
        EXPECT_NE(hbonds[0].donor->name.find("N"), std::string::npos);
      }

      // acceptor chain A, res 46, atom name contains "O"
      if (hbonds[0].acceptor) {
        EXPECT_EQ(atom_chain[hbonds[0].acceptor], "A");
        EXPECT_EQ(atom_resno[hbonds[0].acceptor], 46);
        EXPECT_NE(hbonds[0].acceptor->name.find("O"), std::string::npos);
      }
    }

    if (hbonds.size() >= 2) {
      EXPECT_NEAR(hbonds[1].dist, 2.83621, 1e-4);
      EXPECT_NEAR(hbonds[1].angle_1, 93.8276, 1e-4);
      EXPECT_NEAR(hbonds[1].angle_2, -1.0, 1e-4);
      EXPECT_NEAR(hbonds[1].angle_3, 158.257, 0.001);
      EXPECT_EQ(hbonds[1].ligand_atom_is_donor, false);
      EXPECT_EQ(hbonds[1].bond_has_hydrogen_flag, true);

      // donor chain A, res 47, atom name contains "N"
      if (hbonds[1].donor) {
        EXPECT_EQ(atom_chain[hbonds[1].donor], "A");
        EXPECT_EQ(atom_resno[hbonds[1].donor], 47);
        EXPECT_NE(hbonds[1].donor->name.find("N"), std::string::npos);
      }

      // acceptor chain A, res 46, atom name contains "OG"
      if (hbonds[1].acceptor) {
        EXPECT_EQ(atom_chain[hbonds[1].acceptor], "A");
        EXPECT_EQ(atom_resno[hbonds[1].acceptor], 46);
        EXPECT_NE(hbonds[1].acceptor->name.find("OG"), std::string::npos);
      }
    }
  }

  // Case 2: Use a non-polar residue (ALA 93) — fewer H-bonds expected
  {
    auto sel_all = select_atoms(model, "*", 0, "", 99999, "", true);
    auto sel_lig = select_atoms(model, "A", 93, "", 93, "");

    double max_dist = 3.8;

    auto hbonds = coot::get_mcdonald_and_thornton_gemmi(
        sel_lig, sel_all, model, geom, 0, max_dist);

    EXPECT_EQ(hbonds.size(), 1u);
  }

  // Case 3: Very tight max_dist — should find nothing
  {
    auto sel_all = select_atoms(model, "*", 0, "", 99999, "", true);
    auto sel_lig = select_atoms(model, "A", 46, "", 46, "");

    double max_dist = 0.1;

    auto hbonds = coot::get_mcdonald_and_thornton_gemmi(
        sel_lig, sel_all, model, geom, 0, max_dist);

    EXPECT_EQ(hbonds.size(), 0u);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
