#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, b_factor_histogram) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model& model = st.models[0];

    // Case 1: Select all atoms in chain A
    {
        coot::b_factor_histogram hist =
            coot::b_factor_histogram::b_factor_histogram_gemmi(model, 'A');

        EXPECT_EQ(hist.n_atoms, 2106);
        EXPECT_NEAR(hist.b_max, 105.07f, 1e-4);
        EXPECT_EQ(hist.n_bins, 40);
        EXPECT_EQ(hist.b_vector.size(), 40u);

        // Bin counts (first 5 bins)
        EXPECT_EQ(hist.b_vector[0].size(), 0u);
        EXPECT_EQ(hist.b_vector[1].size(), 0u);
        EXPECT_EQ(hist.b_vector[2].size(), 0u);
        EXPECT_EQ(hist.b_vector[3].size(), 0u);
        EXPECT_EQ(hist.b_vector[4].size(), 0u);
    }

    // Case 2: Empty selection (chain that doesn't exist)
    {
        coot::b_factor_histogram hist =
            coot::b_factor_histogram::b_factor_histogram_gemmi(model, 'Z');

        EXPECT_EQ(hist.n_atoms, 0);
        // NOTE: n_bins and b_max are uninitialized when n_atoms==0 (early return in constructor).
        // The oracle printed 0 and -1 respectively, but these are UB — not asserted.
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}