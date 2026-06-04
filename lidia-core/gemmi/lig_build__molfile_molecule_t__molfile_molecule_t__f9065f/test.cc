#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

TEST(OracleTest, lig_build_molfile_molecule_t) {
    // Case 1: default construction — the function under test
    {
        lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t::molfile_molecule_t_gemmi();
        EXPECT_EQ(mol.atoms.size(), 0u);
        EXPECT_EQ(mol.bonds.size(), 0u);
    }

    // Case 2: verify state after adding an atom (uses same default ctor + add_atom)
    {
        lig_build::molfile_molecule_t mol = lig_build::molfile_molecule_t::molfile_molecule_t_gemmi();
        gemmi::Position pos(1.0, 2.0, 3.0);
        lig_build::molfile_atom_t atom(pos, "C", " C  ");
        EXPECT_EQ(mol.atoms.size(), 0u);
        mol.add_molfile_atom(atom);
        EXPECT_EQ(mol.atoms.size(), 1u);
        EXPECT_EQ(mol.bonds.size(), 0u);
    }

    // Case 3: empty molecule from default ctor — confirm no side effects
    {
        lig_build::molfile_molecule_t mol1 = lig_build::molfile_molecule_t::molfile_molecule_t_gemmi();
        lig_build::molfile_molecule_t mol2 = lig_build::molfile_molecule_t::molfile_molecule_t_gemmi();
        EXPECT_EQ(mol1.atoms.size(), 0u);
        EXPECT_EQ(mol2.atoms.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}