#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, install_ligand) {
    // Case 1: Install a ligand with one fragment, one residue, two atoms
    {
        coot::ligand lig;

        coot::minimol::molecule mol;
        coot::minimol::fragment frag("A");
        coot::minimol::residue res(1, "GLY");
        res.addatom("N", "N", 10.0f, 20.0f, 30.0f, "", 1.0f, 20.0f);
        res.addatom("CA", "C", 11.0f, 21.0f, 31.0f, "", 1.0f, 20.0f);
        frag.residues.push_back(res);
        mol.fragments.push_back(frag);

        EXPECT_EQ(res.n_atoms(), 2u);
        EXPECT_EQ(res.seqnum, 1);
        EXPECT_EQ(res.name, "GLY");
        EXPECT_EQ(frag.fragment_id, "A");

        EXPECT_EQ(lig.initial_ligand.size(), 0);

        lig.install_ligand_gemmi(mol);

        EXPECT_EQ(lig.initial_ligand.size(), 1);
    }

    // Case 2: Install a second ligand — verifies multiple installments work
    {
        coot::ligand lig;

        // First install (from Case 1)
        {
            coot::minimol::molecule mol;
            coot::minimol::fragment frag("A");
            coot::minimol::residue res(1, "GLY");
            res.addatom("N", "N", 10.0f, 20.0f, 30.0f, "", 1.0f, 20.0f);
            res.addatom("CA", "C", 11.0f, 21.0f, 31.0f, "", 1.0f, 20.0f);
            frag.residues.push_back(res);
            mol.fragments.push_back(frag);
            lig.install_ligand_gemmi(mol);
        }

        coot::minimol::molecule mol2;
        coot::minimol::fragment frag2("B");
        coot::minimol::residue res2(42, "ALA");
        res2.addatom("CB", "C", 5.0f, 6.0f, 7.0f, "", 0.5f, 10.0f);
        frag2.residues.push_back(res2);
        mol2.fragments.push_back(frag2);

        EXPECT_EQ(res2.n_atoms(), 1u);
        EXPECT_EQ(res2.seqnum, 42);
        EXPECT_EQ(res2.name, "ALA");
        EXPECT_EQ(frag2.fragment_id, "B");

        EXPECT_EQ(lig.initial_ligand.size(), 1);

        lig.install_ligand_gemmi(mol2);

        EXPECT_EQ(lig.initial_ligand.size(), 2);
    }

    // Case 3: Install an empty molecule (edge case)
    {
        coot::ligand lig;

        // First two installs
        {
            coot::minimol::molecule mol;
            coot::minimol::fragment frag("A");
            coot::minimol::residue res(1, "GLY");
            res.addatom("N", "N", 10.0f, 20.0f, 30.0f, "", 1.0f, 20.0f);
            res.addatom("CA", "C", 11.0f, 21.0f, 31.0f, "", 1.0f, 20.0f);
            frag.residues.push_back(res);
            mol.fragments.push_back(frag);
            lig.install_ligand_gemmi(mol);
        }
        {
            coot::minimol::molecule mol2;
            coot::minimol::fragment frag2("B");
            coot::minimol::residue res2(42, "ALA");
            res2.addatom("CB", "C", 5.0f, 6.0f, 7.0f, "", 0.5f, 10.0f);
            frag2.residues.push_back(res2);
            mol2.fragments.push_back(frag2);
            lig.install_ligand_gemmi(mol2);
        }

        coot::minimol::molecule mol_empty;
        // No fragments, residues, or atoms added

        EXPECT_EQ(mol_empty.fragments.size(), 0u);

        EXPECT_EQ(lig.initial_ligand.size(), 2);

        lig.install_ligand_gemmi(mol_empty);

        EXPECT_EQ(lig.initial_ligand.size(), 3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
