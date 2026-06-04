#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string pad_to_4(const std::string& s) {
    return s + std::string(4 - s.size(), ' ');
}

TEST(OracleTest, chem_feat_clust) {
    // Load the PDB directly with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Case 1: Valid ligand (protein residue) with a HOH water nearby
    {
        coot::residue_spec_t lig_res("A", 10, "");
        std::vector<coot::residue_spec_t> water_specs;
        water_specs.push_back(coot::residue_spec_t("B", 1, ""));

        coot::chem_feat_clust::chem_feat_solvated_ligand_spec_gemmi lig_spec(
            lig_res, water_specs, &st, 0);

        std::vector<coot::residue_spec_t> protein_residues;
        std::vector<coot::chem_feat_clust::chem_feat_solvated_ligand_spec_gemmi> ligand_specs;
        ligand_specs.push_back(lig_spec);

        double water_dist_cutoff = 1000.0;

        coot::chem_feat_clust::chem_feat_clust_gemmi cfc(
            protein_residues, ligand_specs, water_dist_cutoff, true);

        // Call fill_waters to populate water_positions
        coot::chem_feat_clust::fill_waters_gemmi(cfc);

        EXPECT_TRUE(cfc.setup_success);
        EXPECT_EQ(cfc.water_positions.size(), 1u);

        EXPECT_EQ(cfc.water_positions[0].imol, 0);
        EXPECT_EQ(cfc.water_positions[0].ligand_idx, 0u);
        EXPECT_EQ(cfc.water_positions[0].water_spec_idx, 0u);
        EXPECT_NE(cfc.water_positions[0].atom_p, nullptr);
        EXPECT_EQ(pad_to_4(cfc.water_positions[0].atom_p->padded_name()), " O  ");
    }

    // Case 2: No waters specified — water_positions should remain empty
    {
        coot::residue_spec_t lig_res("A", 10, "");
        std::vector<coot::residue_spec_t> water_specs; // empty

        coot::chem_feat_clust::chem_feat_solvated_ligand_spec_gemmi lig_spec(
            lig_res, water_specs, &st, 0);

        std::vector<coot::residue_spec_t> protein_residues;
        std::vector<coot::chem_feat_clust::chem_feat_solvated_ligand_spec_gemmi> ligand_specs;
        ligand_specs.push_back(lig_spec);

        coot::chem_feat_clust::chem_feat_clust_gemmi cfc(
            protein_residues, ligand_specs, 4.2, true);

        // Call fill_waters to populate water_positions
        coot::chem_feat_clust::fill_waters_gemmi(cfc);

        EXPECT_TRUE(cfc.setup_success);
        EXPECT_EQ(cfc.water_positions.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}