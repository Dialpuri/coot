#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find residue by chain name and sequence number
static gemmi::Residue* find_residue(gemmi::Model& model,
                                     const std::string& chain_name,
                                     int seq_num) {
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num)
                    return &res;
            }
        }
    }
    return nullptr;
}

TEST(OracleTestGemmi, normalize_density_blocks_all_positive) {
    // Read structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    // Get residue //A/10
    gemmi::Residue* res = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res, nullptr);

    // === Case 1: All positive values (original oracle case) ===
    {
        std::map<gemmi::Residue*, coot::density_box_t_gemmi> density_block_map_cache;
        double mn_scale_for_normalized_density = 1.0;

        // n_steps=1 => box size = (2*1+1)^3 = 27
        int n_steps = 1;
        int nnn = (2 * n_steps + 1) * (2 * n_steps + 1) * (2 * n_steps + 1);

        ASSERT_EQ(nnn, 27);

        float* box_data = new float[nnn];
        // Fill with positive values: all 10.0
        for (int i = 0; i < nnn; i++) {
            box_data[i] = 10.0f;
        }

        // Create density_box_t_gemmi from raw data
        coot::density_box_t_gemmi db(box_data, res, n_steps);

        // Compute sum before normalization
        double sum_before = 0;
        int count_before = 0;
        for (int i = 0; i < nnn; i++) {
            if (box_data[i] > 0.0) {
                sum_before += box_data[i];
                count_before++;
            }
        }
        double avg_before = (count_before > 0) ? sum_before / count_before : 0.0;

        EXPECT_NEAR(sum_before, 270, 1e-4);
        EXPECT_NEAR(avg_before, 10, 1e-4);

        // Populate the map
        density_block_map_cache[res] = db;

        // Call the function
        coot::side_chain_densities::normalize_density_blocks_gemmi(density_block_map_cache, mn_scale_for_normalized_density);

        // Compute sum after normalization
        double sum_after = 0;
        int count_after = 0;
        coot::density_box_t_gemmi& db_after = density_block_map_cache[res];
        for (int i = 0; i < nnn; i++) {
            if (db_after[i] > 0.0) {
                sum_after += db_after[i];
                count_after++;
            }
        }
        double avg_after = (count_after > 0) ? sum_after / count_after : 0.0;

        EXPECT_NEAR(sum_after, 27, 1e-4);
        EXPECT_NEAR(avg_after, 1, 1e-4);

        // Cleanup
        delete[] box_data;
    }
}

TEST(OracleTestGemmi, normalize_density_blocks_all_non_positive) {
    // Read structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    // Get residue //A/10
    gemmi::Residue* res = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res, nullptr);

    // === Case 2: All non-positive values (n_grid_pts == 0 branch) ===
    {
        std::map<gemmi::Residue*, coot::density_box_t_gemmi> density_block_map_cache;
        double mn_scale_for_normalized_density = 1.0;

        // n_steps=1 => box size = (2*1+1)^3 = 27
        int n_steps = 1;
        int nnn = (2 * n_steps + 1) * (2 * n_steps + 1) * (2 * n_steps + 1);

        float* box_data = new float[nnn];
        // Fill with negative values: all -5.0
        for (int i = 0; i < nnn; i++) {
            box_data[i] = -5.0f;
        }

        // Create density_box_t_gemmi from raw data
        coot::density_box_t_gemmi db(box_data, res, n_steps);

        // Compute sum before normalization (only positive values count)
        double sum_before = 0;
        int count_before = 0;
        for (int i = 0; i < nnn; i++) {
            if (box_data[i] > 0.0) {
                sum_before += box_data[i];
                count_before++;
            }
        }

        // No positive values, so n_grid_pts should be 0
        EXPECT_EQ(count_before, 0);
        EXPECT_NEAR(sum_before, 0, 1e-4);

        // Populate the map
        density_block_map_cache[res] = db;

        // Call the function — should NOT modify data (n_grid_pts == 0)
        coot::side_chain_densities::normalize_density_blocks_gemmi(density_block_map_cache, mn_scale_for_normalized_density);

        // Check that values are unchanged (still -5.0)
        coot::density_box_t_gemmi& db_after = density_block_map_cache[res];
        for (int i = 0; i < nnn; i++) {
            EXPECT_NEAR(db_after[i], -5.0f, 1e-4);
        }

        // Cleanup
        delete[] box_data;
    }
}

TEST(OracleTestGemmi, normalize_density_blocks_mixed_values) {
    // Read structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    // Get residue //A/10
    gemmi::Residue* res = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res, nullptr);

    // === Case 3: Mixed positive and non-positive values ===
    {
        std::map<gemmi::Residue*, coot::density_box_t_gemmi> density_block_map_cache;
        double mn_scale_for_normalized_density = 1.0;

        // n_steps=1 => box size = (2*1+1)^3 = 27
        int n_steps = 1;
        int nnn = (2 * n_steps + 1) * (2 * n_steps + 1) * (2 * n_steps + 1);

        float* box_data = new float[nnn];
        // Fill first half with positive 20.0, second half with -3.0
        for (int i = 0; i < nnn; i++) {
            if (i < 13)
                box_data[i] = 20.0f;
            else
                box_data[i] = -3.0f;
        }

        // Create density_box_t_gemmi from raw data
        coot::density_box_t_gemmi db(box_data, res, n_steps);

        // Before: sum of positives = 13 * 20 = 260, avg of positives = 260/13 = 20
        double sum_before = 0;
        int count_before = 0;
        for (int i = 0; i < nnn; i++) {
            if (box_data[i] > 0.0) {
                sum_before += box_data[i];
                count_before++;
            }
        }
        EXPECT_EQ(count_before, 13);
        EXPECT_NEAR(sum_before, 260, 1e-4);

        // Populate the map
        density_block_map_cache[res] = db;

        // Call the function
        // av = 260/13 = 20, sc = 1.0/20 = 0.05
        coot::side_chain_densities::normalize_density_blocks_gemmi(density_block_map_cache, mn_scale_for_normalized_density);

        // After: positive values scaled by 0.05 => 20 * 0.05 = 1.0 (for first 13)
        // negative values scaled by 0.05 => -3.0 * 0.05 = -0.15 (for last 14)
        // sum after positives: 13 * 1.0 = 13
        // avg after positives: 13/13 = 1.0
        coot::density_box_t_gemmi& db_after = density_block_map_cache[res];
        for (int i = 0; i < 13; i++) {
            EXPECT_NEAR(db_after[i], 1.0, 1e-4);
        }
        for (int i = 13; i < nnn; i++) {
            EXPECT_NEAR(db_after[i], -0.15, 1e-4);
        }

        double sum_after = 0;
        int count_after = 0;
        for (int i = 0; i < nnn; i++) {
            if (db_after[i] > 0.0) {
                sum_after += db_after[i];
                count_after++;
            }
        }
        double avg_after = (count_after > 0) ? sum_after / count_after : 0.0;

        EXPECT_NEAR(sum_after, 13, 1e-4);
        EXPECT_NEAR(avg_after, 1, 1e-4);

        // Cleanup
        delete[] box_data;
    }
}

TEST(OracleTestGemmi, normalize_density_blocks_custom_scale) {
    // Read structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);

    // Get residue //A/10
    gemmi::Residue* res = find_residue(st.models[0], "A", 10);
    ASSERT_NE(res, nullptr);

    // === Case 4: Custom scale factor ===
    {
        std::map<gemmi::Residue*, coot::density_box_t_gemmi> density_block_map_cache;
        double mn_scale_for_normalized_density = 2.0; // custom scale

        int n_steps = 1;
        int nnn = (2 * n_steps + 1) * (2 * n_steps + 1) * (2 * n_steps + 1);

        float* box_data = new float[nnn];
        for (int i = 0; i < nnn; i++) {
            box_data[i] = 10.0f;
        }

        coot::density_box_t_gemmi db(box_data, res, n_steps);

        density_block_map_cache[res] = db;

        // av = 270/27 = 10, sc = 2.0/10 = 0.2
        coot::side_chain_densities::normalize_density_blocks_gemmi(density_block_map_cache, mn_scale_for_normalized_density);

        coot::density_box_t_gemmi& db_after = density_block_map_cache[res];
        for (int i = 0; i < nnn; i++) {
            EXPECT_NEAR(db_after[i], 2.0, 1e-4); // 10 * 0.2 = 2.0
        }

        double sum_after = 0;
        int count_after = 0;
        for (int i = 0; i < nnn; i++) {
            if (db_after[i] > 0.0) {
                sum_after += db_after[i];
                count_after++;
            }
        }
        EXPECT_NEAR(sum_after, 54, 1e-4); // 27 * 2 = 54

        delete[] box_data;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}