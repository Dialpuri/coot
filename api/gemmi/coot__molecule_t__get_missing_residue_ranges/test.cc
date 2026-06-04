#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_missing_residue_ranges) {
    {
        // case: example.pdb — the original oracle case
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::vector<coot::residue_range_t> ranges = coot::get_missing_residue_ranges_gemmi(st.models[0]);

        EXPECT_EQ(ranges.size(), 3u);

        EXPECT_EQ(ranges[0].chain_id, "");
        EXPECT_EQ(ranges[0].res_no_start, 37);
        EXPECT_EQ(ranges[0].res_no_end, 45);

        EXPECT_EQ(ranges[1].chain_id, "");
        EXPECT_EQ(ranges[1].res_no_start, 73);
        EXPECT_EQ(ranges[1].res_no_end, 75);

        EXPECT_EQ(ranges[2].chain_id, "");
        EXPECT_EQ(ranges[2].res_no_start, 147);
        EXPECT_EQ(ranges[2].res_no_end, 165);
    }
    {
        // case: empty molecule (no model loaded) — exercises null-pointer branches
        gemmi::Model empty_model;

        std::vector<coot::residue_range_t> ranges = coot::get_missing_residue_ranges_gemmi(empty_model);

        EXPECT_TRUE(ranges.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}