#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_covalent_bonds_by_distance) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Case 1: valid residue with neighbors (residue 10 in chain A — should find 2 peptide bonds)
    {
        coot::residue_spec_t ligand_spec("A", 10, "");

        std::vector<pli::fle_ligand_bond_t> bonds =
            pli::get_covalent_bonds_by_distance_gemmi(model, ligand_spec);

        EXPECT_EQ(bonds.size(), 2u);

        // bond_0: lig_atom=A/10/ N  env_atom=A/9/ C  dist_sq=1.78285
        EXPECT_EQ(bonds[0].ligand_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[0].ligand_atom_spec.res_no, 10);
        EXPECT_EQ(bonds[0].ligand_atom_spec.atom_name, " N  ");
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.res_no, 9);
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.atom_name, " C  ");
        EXPECT_NEAR(bonds[0].bond_length, 1.78285, 1e-4);

        // bond_1: lig_atom=A/10/ C  env_atom=A/11/ N  dist_sq=1.79079
        EXPECT_EQ(bonds[1].ligand_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[1].ligand_atom_spec.res_no, 10);
        EXPECT_EQ(bonds[1].ligand_atom_spec.atom_name, " C  ");
        EXPECT_EQ(bonds[1].interacting_residue_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[1].interacting_residue_atom_spec.res_no, 11);
        EXPECT_EQ(bonds[1].interacting_residue_atom_spec.atom_name, " N  ");
        EXPECT_NEAR(bonds[1].bond_length, 1.79079, 1e-4);
    }

    // Case 2: terminal residue (residue 1 in chain A — should find 1 peptide bond)
    {
        coot::residue_spec_t ligand_spec2("A", 1, "");

        std::vector<pli::fle_ligand_bond_t> bonds =
            pli::get_covalent_bonds_by_distance_gemmi(model, ligand_spec2);

        EXPECT_EQ(bonds.size(), 1u);

        // bond_0: lig_atom=A/1/ C  env_atom=A/2/ N  dist_sq=1.79907
        EXPECT_EQ(bonds[0].ligand_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[0].ligand_atom_spec.res_no, 1);
        EXPECT_EQ(bonds[0].ligand_atom_spec.atom_name, " C  ");
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.chain_id, "A");
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.res_no, 2);
        EXPECT_EQ(bonds[0].interacting_residue_atom_spec.atom_name, " N  ");
        EXPECT_NEAR(bonds[0].bond_length, 1.79907, 1e-4);
    }

    // Case 3: water residue (HOH in chain B, res 1 — should find 0 covalent bonds)
    {
        coot::residue_spec_t ligand_spec3("B", 1, "");

        std::vector<pli::fle_ligand_bond_t> bonds =
            pli::get_covalent_bonds_by_distance_gemmi(model, ligand_spec3);

        EXPECT_EQ(bonds.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}