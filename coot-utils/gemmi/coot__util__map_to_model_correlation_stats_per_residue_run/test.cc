#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

// We still need molecules_container_t to load the MTZ and get the xmap
#include "api/molecules-container.hh"

TEST(OracleTest, map_to_model_correlation_stats_per_residue_run) {
    // Load MTZ to get an electron density map (xmap is clipper-based, so we use molecules_container)
    molecules_container_t mc;
    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    EXPECT_GE(imol, 0);

    // Load MTZ to get an electron density map into the molecule
    mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz", "FWT", "PHWT", "", false, false);

    // Find the molecule with the xmap
    int imol_map = -1;
    int n_mol = mc.get_number_of_molecules();
    for (int i = 0; i < n_mol; i++) {
        if (!mc[i].xmap.is_null()) {
            imol_map = i;
            break;
        }
    }
    ASSERT_GE(imol_map, 0);

    clipper::Xmap<float> xmap = mc[imol_map].xmap;

    // Load the same structure with gemmi for the model
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());
    const gemmi::Model &model = st.models[0];

    std::string chain_id = "A";
    unsigned int n_residues_per_blob = 3;
    bool exclude_CON = true;
    float atom_mask_radius = 2.8f;
    float NOC_mask_radius = 1.8f;

    // Normal case: valid chain
    {
        auto result = coot::util::map_to_model_correlation_stats_per_residue_run_gemmi(
            model, chain_id, xmap, n_residues_per_blob, exclude_CON, atom_mask_radius, NOC_mask_radius);

        EXPECT_EQ(result.first.size(), 273u);
        EXPECT_EQ(result.second.size(), 273u);

        ASSERT_FALSE(result.first.empty());
        auto it_all = result.first.begin();
        EXPECT_EQ(it_all->first.chain_id, "A");
        EXPECT_EQ(it_all->first.res_no, 2);
        EXPECT_EQ(it_all->second.n, 5354);
        if (it_all->second.n > 0) {
            EXPECT_NEAR(it_all->second.correlation(), 0.772009, 1e-4);
        }

        ASSERT_FALSE(result.second.empty());
        auto it_sc = result.second.begin();
        EXPECT_EQ(it_sc->first.chain_id, "A");
        EXPECT_EQ(it_sc->first.res_no, 2);
        EXPECT_EQ(it_sc->second.n, 7086);
        if (it_sc->second.n > 0) {
            EXPECT_NEAR(it_sc->second.correlation(), 0.619276, 1e-4);
        }
    }

    // Edge case: invalid chain ID
    {
        std::string bad_chain = "ZZZ";
        auto result2 = coot::util::map_to_model_correlation_stats_per_residue_run_gemmi(
            model, bad_chain, xmap, n_residues_per_blob, exclude_CON, atom_mask_radius, NOC_mask_radius);
        EXPECT_EQ(result2.first.size(), 0u);
        EXPECT_EQ(result2.second.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}