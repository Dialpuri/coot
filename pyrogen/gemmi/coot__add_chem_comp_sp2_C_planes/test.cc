#define LIBCOOTAPI_BUILD
#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/RDGeneral/Invariant.h>
#include "function.hh"

using namespace RDKit;

TEST(OracleTest, add_chem_comp_sp2_C_planes) {
  // Case 1: valid ketone pattern CC(=O)C matches AC(=O)C
  {
    ROMol *mol_ptr = SmilesToMol("CC(=O)C");
    ASSERT_NE(mol_ptr, nullptr);

    mol_ptr->getAtomWithIdx(0)->setProp("name", "C1");
    mol_ptr->getAtomWithIdx(1)->setProp("name", "C2");
    mol_ptr->getAtomWithIdx(2)->setProp("name", "O1");
    mol_ptr->getAtomWithIdx(3)->setProp("name", "C3");

    const ROMol &mol = *mol_ptr;
    EXPECT_EQ(mol.getNumAtoms(), 4);

    coot::dictionary_residue_restraints_t restraints;
    int planes_before = restraints.plane_restraint.size();

    coot::add_chem_comp_sp2_C_planes_gemmi(mol, &restraints);

    int planes_after = restraints.plane_restraint.size();

    EXPECT_EQ(planes_before, 0);
    EXPECT_EQ(planes_after, 1);
    EXPECT_EQ(restraints.plane_restraint[0].plane_id, "plane-sp2-N-1");

    delete mol_ptr;
  }

  // Case 2: no ketone — plain alkane, should produce no planes
  {
    ROMol *mol_ptr = SmilesToMol("CCCC");
    ASSERT_NE(mol_ptr, nullptr);

    mol_ptr->getAtomWithIdx(0)->setProp("name", "C1");
    mol_ptr->getAtomWithIdx(1)->setProp("name", "C2");
    mol_ptr->getAtomWithIdx(2)->setProp("name", "C3");
    mol_ptr->getAtomWithIdx(3)->setProp("name", "C4");

    const ROMol &mol = *mol_ptr;
    EXPECT_EQ(mol.getNumAtoms(), 4);

    coot::dictionary_residue_restraints_t restraints;
    int planes_before = restraints.plane_restraint.size();

    coot::add_chem_comp_sp2_C_planes_gemmi(mol, &restraints);

    int planes_after = restraints.plane_restraint.size();

    EXPECT_EQ(planes_before, 0);
    EXPECT_EQ(planes_after, 0);

    delete mol_ptr;
  }

  // Case 3: acetic acid CC(=O)O — should still match the pattern
  {
    ROMol *mol_ptr = SmilesToMol("CC(=O)O");
    ASSERT_NE(mol_ptr, nullptr);

    mol_ptr->getAtomWithIdx(0)->setProp("name", "C1");
    mol_ptr->getAtomWithIdx(1)->setProp("name", "C2");
    mol_ptr->getAtomWithIdx(2)->setProp("name", "O1");
    mol_ptr->getAtomWithIdx(3)->setProp("name", "OH");

    const ROMol &mol = *mol_ptr;
    EXPECT_EQ(mol.getNumAtoms(), 4);

    coot::dictionary_residue_restraints_t restraints;
    int planes_before = restraints.plane_restraint.size();

    coot::add_chem_comp_sp2_C_planes_gemmi(mol, &restraints);

    int planes_after = restraints.plane_restraint.size();

    EXPECT_EQ(planes_before, 0);
    EXPECT_EQ(planes_after, 1);

    delete mol_ptr;
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}