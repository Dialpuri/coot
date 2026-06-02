#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_cablam_like_geometry_stats) {
    // Case: normal molecule
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<coot::cablam_like_geometry_stats_t_gemmi> stats = coot::get_cablam_like_geometry_stats_gemmi(st);

        EXPECT_EQ(stats.size(), 258u);

        // entry 0
        EXPECT_NE(stats[0].residue, nullptr);
        std::string residue_0 = std::string("//") +
            stats[0].chain_id + "/" +
            std::to_string(stats[0].residue->seqid.num.value);
        EXPECT_EQ(residue_0, "//A/2");
        EXPECT_NEAR(stats[0].dp_prev_to_mid, 1.43555, 1e-4);
        EXPECT_NEAR(stats[0].dp_next_to_mid, 0.66697, 1e-4);
        EXPECT_NEAR(stats[0].dist_proj_point_prev_to_next, 4.27661, 1e-4);

        // entry 1
        EXPECT_NE(stats[1].residue, nullptr);
        std::string residue_1 = std::string("//") +
            stats[1].chain_id + "/" +
            std::to_string(stats[1].residue->seqid.num.value);
        EXPECT_EQ(residue_1, "//A/3");
        EXPECT_NEAR(stats[1].dp_prev_to_mid, 0.66697, 1e-4);
        EXPECT_NEAR(stats[1].dp_next_to_mid, -3.02665, 1e-4);
        EXPECT_NEAR(stats[1].dist_proj_point_prev_to_next, 5.3575, 1e-4);

        // entry 2
        EXPECT_NE(stats[2].residue, nullptr);
        std::string residue_2 = std::string("//") +
            stats[2].chain_id + "/" +
            std::to_string(stats[2].residue->seqid.num.value);
        EXPECT_EQ(residue_2, "//A/4");
        EXPECT_NEAR(stats[2].dp_prev_to_mid, -3.02665, 1e-4);
        EXPECT_NEAR(stats[2].dp_next_to_mid, -3.00725, 1e-4);
        EXPECT_NEAR(stats[2].dist_proj_point_prev_to_next, 6.64848, 1e-4);
    }

    // Case: empty molecule (no models)
    {
        gemmi::Structure st_empty;
        std::vector<coot::cablam_like_geometry_stats_t_gemmi> stats_empty = coot::get_cablam_like_geometry_stats_gemmi(st_empty);
        EXPECT_EQ(stats_empty.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}