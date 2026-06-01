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

TEST(OracleTest, add_chem_comp_deloc_planes) {
    // Case 1: Create a molecule with a carboxylate group (matches pattern *C(=O)[O;H])
    // Using acetic acid: CC(=O)O  — the C(=O)O group should match
    {
        ROMol *mol1 = SmilesToMol("CC(=O)O");
        ASSERT_NE(mol1, nullptr);

        // Set "name" property on each atom (required by add_chem_comp_deloc_planes)
        mol1->getAtomWithIdx(0)->setProp("name", "C1");
        mol1->getAtomWithIdx(1)->setProp("name", "C2");
        mol1->getAtomWithIdx(2)->setProp("name", "O1");
        mol1->getAtomWithIdx(3)->setProp("name", "OH");

        coot::dictionary_residue_restraints_t restraints1;

        EXPECT_EQ(restraints1.plane_restraint.size(), 0u);

        coot::add_chem_comp_deloc_planes_gemmi(*mol1, &restraints1);

        EXPECT_EQ(restraints1.plane_restraint.size(), 2u);

        {
            const coot::dict_plane_restraint_t &pr = restraints1.plane_restraint[0];
            EXPECT_EQ(pr.plane_id, "plane-deloc-1");
            EXPECT_EQ(pr.n_atoms(), 4);
            EXPECT_EQ(pr.atom_id(0), " C1 ");
            EXPECT_EQ(pr.atom_id(1), " C2 ");
            EXPECT_EQ(pr.atom_id(2), " O1 ");
            EXPECT_EQ(pr.atom_id(3), " OH ");
        }
        {
            const coot::dict_plane_restraint_t &pr = restraints1.plane_restraint[1];
            EXPECT_EQ(pr.plane_id, "plane-deloc-2");
            EXPECT_EQ(pr.n_atoms(), 4);
            EXPECT_EQ(pr.atom_id(0), " C1 ");
            EXPECT_EQ(pr.atom_id(1), " C2 ");
            EXPECT_EQ(pr.atom_id(2), " OH ");
            EXPECT_EQ(pr.atom_id(3), " O1 ");
        }
        delete mol1;
    }

    // Case 2 and 3 from oracle crashed (SmilesToMol("C") has only 1 explicit atom).
    // Only Case 1 produced observed output.
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}