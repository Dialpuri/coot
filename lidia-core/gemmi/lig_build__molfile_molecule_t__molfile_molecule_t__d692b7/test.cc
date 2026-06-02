#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Case 1: Valid restraints with atoms and bonds
TEST(MolfileMoleculeGemmiTest, ValidRestraints) {
    coot::dictionary_residue_restraints_t restraints;

    // Add two atoms: CA and CB — use 4-char format matching bond_restraint constructor output
    coot::dict_atom atom_ca;
    atom_ca.atom_id = " CA ";
    atom_ca.type_symbol = "C";
    restraints.atom_info.push_back(atom_ca);

    coot::dict_atom atom_cb;
    atom_cb.atom_id = " CB ";
    atom_cb.type_symbol = "C";
    restraints.atom_info.push_back(atom_cb);

    // Add one bond: CA-CB (single) — constructor pads " CA " → " CA " (substr 0,4)
    coot::dict_bond_restraint_t bond_ca_cb(" CA ", " CB ", "single");
    restraints.bond_restraint.push_back(bond_ca_cb);

    EXPECT_EQ(restraints.atom_info.size(), 2u);
    EXPECT_EQ(restraints.bond_restraint.size(), 1u);
    EXPECT_EQ(restraints.atom_info[0].atom_id, " CA ");
    EXPECT_EQ(restraints.atom_info[1].atom_id, " CB ");

    lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t_gemmi(restraints);

    EXPECT_EQ(mol.atoms.size(), 3u);  // blank + CA + CB
    EXPECT_EQ(mol.bonds.size(), 1u);
    EXPECT_EQ(mol.atoms[0].name, "");
    EXPECT_EQ(mol.atoms[0].element, "");
    EXPECT_EQ(mol.atoms[1].name, " CA ");
    EXPECT_EQ(mol.atoms[1].element, "C");
    EXPECT_EQ(mol.atoms[2].name, " CB ");
    EXPECT_EQ(mol.atoms[2].element, "C");
    EXPECT_EQ(mol.bonds[0].index_1, 1);
    EXPECT_EQ(mol.bonds[0].index_2, 2);
}

// Case 2: Empty restraints — no atoms, no bonds
TEST(MolfileMoleculeGemmiTest, EmptyRestraints) {
    coot::dictionary_residue_restraints_t restraints;

    EXPECT_EQ(restraints.atom_info.size(), 0u);
    EXPECT_EQ(restraints.bond_restraint.size(), 0u);

    lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t_gemmi(restraints);

    EXPECT_EQ(mol.atoms.size(), 1u);  // just the blank atom
    EXPECT_EQ(mol.bonds.size(), 0u);
}

// Case 3: Atoms exist but bond references missing atom — bond should be skipped
TEST(MolfileMoleculeGemmiTest, MissingBondAtom) {
    coot::dictionary_residue_restraints_t restraints;

    coot::dict_atom atom_ca;
    atom_ca.atom_id = " CA ";
    atom_ca.type_symbol = "C";
    restraints.atom_info.push_back(atom_ca);

    // Bond references CA and CB, but CB not in atom_info
    coot::dict_bond_restraint_t bond(" CA ", " CB ", "single");
    restraints.bond_restraint.push_back(bond);

    EXPECT_EQ(restraints.atom_info.size(), 1u);
    EXPECT_EQ(restraints.bond_restraint.size(), 1u);
    EXPECT_EQ(restraints.atom_info[0].atom_id, " CA ");

    lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t_gemmi(restraints);

    EXPECT_EQ(mol.atoms.size(), 2u);  // blank + CA
    EXPECT_EQ(mol.bonds.size(), 0u);  // bond skipped because CB missing
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}