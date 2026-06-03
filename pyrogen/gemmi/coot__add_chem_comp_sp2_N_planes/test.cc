#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/MolOps.h>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <utility>
#include <vector>
#include <string>

using namespace RDKit;

TEST(OracleTest, add_chem_comp_sp2_N_planes) {

   // === Oracle Case 1: 6-aminopurine (adenine base) ===
   {
      ROMol *mol = SmilesToMol("Nc1ncnc2ncnn12");
      ASSERT_NE(mol, static_cast<ROMol*>(nullptr));
      MolOps::addHs(*mol);

      std::vector<std::string> names = {"N6", "C6", "N1", "C2", "N3", "C4",
                                         "C5", "N7", "C8", "N9", "N2"};
      for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
         if (i < names.size()) {
            mol->getAtomWithIdx(i)->setProp("name", names[i]);
         } else {
            std::ostringstream oss;
            oss << "H" << i;
            mol->getAtomWithIdx(i)->setProp("name", oss.str());
         }
      }

      coot::dictionary_residue_restraints_t restraints;
      int before = restraints.plane_restraint.size();
      coot::add_chem_comp_sp2_N_planes_gemmi(*mol, &restraints);
      int after = restraints.plane_restraint.size();

      EXPECT_EQ(before, 0);
      EXPECT_EQ(after, 0);

      delete mol;
   }

   // === Oracle Case 2: Ethanol — no N at all ===
   {
      ROMol *mol = SmilesToMol("CCO");
      ASSERT_NE(mol, static_cast<ROMol*>(nullptr));
      MolOps::addHs(*mol);

      mol->getAtomWithIdx(0)->setProp("name", "C1");
      mol->getAtomWithIdx(1)->setProp("name", "C2");
      mol->getAtomWithIdx(2)->setProp("name", "O1");
      for (unsigned int i = 3; i < mol->getNumAtoms(); i++) {
         mol->getAtomWithIdx(i)->setProp("name", "H");
      }

      coot::dictionary_residue_restraints_t restraints;
      int before = restraints.plane_restraint.size();
      coot::add_chem_comp_sp2_N_planes_gemmi(*mol, &restraints);
      int after = restraints.plane_restraint.size();

      EXPECT_EQ(before, 0);
      EXPECT_EQ(after, 0);

      delete mol;
   }

   // === Oracle Case 3: Methylamine — aliphatic NH2, no sp2 C ===
   {
      ROMol *mol = SmilesToMol("CN");
      ASSERT_NE(mol, static_cast<ROMol*>(nullptr));
      MolOps::addHs(*mol);

      mol->getAtomWithIdx(0)->setProp("name", "C1");
      mol->getAtomWithIdx(1)->setProp("name", "N1");
      for (unsigned int i = 2; i < mol->getNumAtoms(); i++) {
         mol->getAtomWithIdx(i)->setProp("name", "H");
      }

      coot::dictionary_residue_restraints_t restraints;
      int before = restraints.plane_restraint.size();
      coot::add_chem_comp_sp2_N_planes_gemmi(*mol, &restraints);
      int after = restraints.plane_restraint.size();

      EXPECT_EQ(before, 0);
      EXPECT_EQ(after, 0);

      delete mol;
   }

   // === Complementary Case: allylamine (C=CN) — sp2 C bonded to NH2 ===
   // This tests a non-aromatic sp2 C (=C) bonded to NH2.
   // The oracle's SMARTS [c,C][N^3;H2]([H])[H] uses N^3, which may
   // have subtle issues with RDKit's valence counting. The oracle
   // observed 0 restraints for all test cases.
   {
      ROMol *mol = SmilesToMol("C=CN");
      ASSERT_NE(mol, static_cast<ROMol*>(nullptr));
      MolOps::addHs(*mol);

      mol->getAtomWithIdx(0)->setProp("name", "C1");
      mol->getAtomWithIdx(1)->setProp("name", "C2");
      mol->getAtomWithIdx(2)->setProp("name", "N1");
      for (unsigned int i = 3; i < mol->getNumAtoms(); i++) {
         std::ostringstream oss;
         oss << "H" << i;
         mol->getAtomWithIdx(i)->setProp("name", oss.str());
      }

      coot::dictionary_residue_restraints_t restraints;
      int before = restraints.plane_restraint.size();
      coot::add_chem_comp_sp2_N_planes_gemmi(*mol, &restraints);
      int after = restraints.plane_restraint.size();

      EXPECT_EQ(before, 0);
      // Consistent with oracle: the N^3 SMARTS pattern doesn't match
      // even when an sp2 C is bonded to NH2.
      EXPECT_EQ(after, 0);

      delete mol;
   }

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}