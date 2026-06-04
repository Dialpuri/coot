#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"

#include <vector>
#include <string>
#include <utility>

TEST(OracleTest, partition_map_by_chain) {

    // Case 1: Original oracle case — protein-ligand model with one chain (chain A)
    {
        molecules_container_t mc;
        mc.geometry_init_standard();

        // Load an MTZ map molecule
        int imol_map = mc.read_mtz(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
            "FWT", "PHWT", "", false, false);
        EXPECT_GE(imol_map, 0);
        EXPECT_EQ(imol_map, 0);

        // Get the xmap from the map molecule
        const clipper::Xmap<float> &xmap = mc[imol_map].xmap;

        EXPECT_FALSE(xmap.is_null());

        EXPECT_NEAR(xmap.cell().a(), 53.512, 1e-4);
        EXPECT_NEAR(xmap.cell().b(), 71.311, 1e-4);
        EXPECT_NEAR(xmap.cell().c(), 72.101, 1e-4);

        // Load structure with gemmi instead of MMDB
        gemmi::Structure st = gemmi::read_structure_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

        EXPECT_EQ(st.models.size(), 1u);
        const gemmi::Model* model_p = &st.models[0];

        // Count chains in the model
        int total_chains = 0;
        for (const auto& model : st.models) {
            total_chains += model.chains.size();
        }
        EXPECT_EQ(total_chains, 1);

        // Call the function
        std::string state_string;
        std::vector<std::pair<std::string, clipper::Xmap<float>>> result =
            coot::util::partition_map_by_chain_gemmi(xmap, model_p, &state_string);

        EXPECT_EQ(result.size(), 1u);
        EXPECT_EQ(state_string, "Constructing map for chain A");

        // Check the partitioned map
        EXPECT_EQ(result[0].first, "A");
        const clipper::Xmap<float> &chain_map = result[0].second;
        EXPECT_FALSE(chain_map.is_null());

        clipper::Grid_sampling gs = chain_map.grid_sampling();
        EXPECT_EQ(gs.nu(), 108);
        EXPECT_EQ(gs.nv(), 144);
        EXPECT_EQ(gs.nw(), 144);
    }

    // Case 2: Complementary case — different PDB file to exercise
    // different code paths (e.g., different chain/backbone atom combinations)
    {
        molecules_container_t mc;
        mc.geometry_init_standard();

        // Load the map
        int imol_map = mc.read_mtz(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
            "FWT", "PHWT", "", false, false);
        EXPECT_GE(imol_map, 0);

        // Get the xmap from the map molecule
        const clipper::Xmap<float> &xmap = mc[imol_map].xmap;

        EXPECT_FALSE(xmap.is_null());

        // Load a different PDB model with gemmi
        gemmi::Structure st = gemmi::read_structure_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        const gemmi::Model* model_p = nullptr;
        if (!st.models.empty()) {
            model_p = &st.models[0];
        }

        // Count chains in the model
        int total_chains = 0;
        for (const auto& model : st.models) {
            total_chains += model.chains.size();
        }

        // Call the function
        std::string state_string;
        std::vector<std::pair<std::string, clipper::Xmap<float>>> result;
        if (model_p) {
            result = coot::util::partition_map_by_chain_gemmi(xmap, model_p, &state_string);
        }

        // Basic checks — result should have at least one entry if there are chains
        if (total_chains > 0) {
            EXPECT_GT(result.size(), 0u);
            EXPECT_FALSE(state_string.empty());

            for (size_t i = 0; i < result.size(); i++) {
                EXPECT_FALSE(result[i].first.empty());
                const clipper::Xmap<float> &chain_map = result[i].second;
                EXPECT_FALSE(chain_map.is_null());
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}