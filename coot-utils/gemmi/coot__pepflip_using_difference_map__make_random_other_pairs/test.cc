#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_random_other_pairs) {
    // Load a PDB into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_EQ(st.models.size(), 1u);

    // Get the first model
    gemmi::Model& model = st.models[0];

    // Build the pepflip_using_difference_map object
    coot::pepflip_using_difference_map pf(&model);

    // Case 1: valid call with n_others=200 (typical)
    {
        int n_others = 200;
        auto result = pf.make_random_other_pairs_gemmi(n_others);
        EXPECT_EQ(result.size(), 193u);
        if (!result.empty()) {
            // Coordinates are random — verify they are finite and non-zero-ish
            EXPECT_FALSE(std::isnan(result[0].first[0]));
            EXPECT_FALSE(std::isnan(result[0].first[1]));
            EXPECT_FALSE(std::isnan(result[0].first[2]));
            EXPECT_FALSE(std::isnan(result[0].second[0]));
            EXPECT_FALSE(std::isnan(result[0].second[1]));
            EXPECT_FALSE(std::isnan(result[0].second[2]));
        }
    }

    // Case 2: n_others=0 should return empty (edge case)
    {
        int n_others = 0;
        auto result = pf.make_random_other_pairs_gemmi(n_others);
        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: very large n_others to get many pairs
    {
        int n_others = 10000;
        auto result = pf.make_random_other_pairs_gemmi(n_others);
        EXPECT_EQ(result.size(), 6218u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}