#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(PepflipTest, no_difference_map_returns_empty) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // case: without difference map (has_valid_xmap == false)
    {
        auto flips = coot::pepflip_using_difference_map::pepflip_using_difference_map_gemmi(
            st, 3.0f, /*has_valid_xmap=*/false);
        EXPECT_EQ(flips.size(), 0);
    }
}

TEST(PepflipTest, empty_structure_returns_empty) {
    gemmi::Structure st; // empty structure
    auto flips = coot::pepflip_using_difference_map::pepflip_using_difference_map_gemmi(
        st, 3.0f, /*has_valid_xmap=*/false);
    EXPECT_EQ(flips.size(), 0);
}

TEST(PepflipTest, valid_coords_non_diff_map_returns_empty) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    // Simulate calling with a non-difference-map molecule: has_valid_xmap == false
    auto flips = coot::pepflip_using_difference_map::pepflip_using_difference_map_gemmi(
        st, 3.0f, /*has_valid_xmap=*/false);
    EXPECT_EQ(flips.size(), 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}