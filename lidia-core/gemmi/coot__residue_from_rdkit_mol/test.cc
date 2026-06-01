#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/Geometry/point.h>

TEST(OracleTest, residue_from_rdkit_mol) {
  using namespace RDKit;

  // Case 1: Valid molecule with named atoms and a manually created conformer
  {
    ROMol *mol = SmilesToMol("O");
    ASSERT_NE(mol, nullptr);

    // Set "name" property on the oxygen atom
    mol->getAtomWithIdx(0)->setProp("name", "O1");

    // Manually create a conformer with coordinates
    unsigned int n_atoms = mol->getNumAtoms();
    Conformer *conf = new Conformer(n_atoms);
    conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
    mol->addConformer(conf);
    int conf_id = conf->getId();

    EXPECT_EQ(n_atoms, 1u);
    EXPECT_EQ(conf_id, 0);

    std::string comp_id = "HOH";
    auto res = coot::residue_from_rdkit_mol_gemmi(*mol, conf_id, comp_id);

    EXPECT_NE(res, nullptr);
    if (res) {
      EXPECT_EQ(res->name, "HOH");
      int natoms = static_cast<int>(res->atoms.size());
      EXPECT_EQ(natoms, 1);
      {
        gemmi::Atom &at = res->atoms[0];
        EXPECT_EQ(at.name, " O1 ");
        EXPECT_EQ(std::string(" ") + at.element.name(), " O");
      }
    }
    delete mol;
  }

  // Case 2: Empty molecule — should return nullptr
  {
    ROMol *mol = new ROMol();
    std::string comp_id = "EMPTY";
    auto res = coot::residue_from_rdkit_mol_gemmi(*mol, 0, comp_id);

    EXPECT_EQ(res, nullptr);
    delete mol;
  }

  // Case 3: Molecule with atoms but one missing "name" property — should skip that atom
  {
    ROMol *mol = SmilesToMol("CC");
    ASSERT_NE(mol, nullptr);
    // Set "name" on only the first atom — second atom will be skipped
    mol->getAtomWithIdx(0)->setProp("name", "C1");

    unsigned int n_atoms = mol->getNumAtoms();
    Conformer *conf = new Conformer(n_atoms);
    conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
    conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
    mol->addConformer(conf);
    int conf_id = conf->getId();

    EXPECT_EQ(n_atoms, 2u);
    EXPECT_EQ(conf_id, 0);

    std::string comp_id = "ETH";
    auto res = coot::residue_from_rdkit_mol_gemmi(*mol, conf_id, comp_id);

    EXPECT_NE(res, nullptr);
    if (res) {
      EXPECT_EQ(res->name, "ETH");
      EXPECT_EQ(res->atoms.size(), 1);
    }
    delete mol;
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}