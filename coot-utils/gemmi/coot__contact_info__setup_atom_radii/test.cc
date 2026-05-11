#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, SetupAtomRadii) {
    auto atom_radii = coot::setup_atom_radii_gemmi();

    // Verify atom_radii size
    EXPECT_EQ(atom_radii.size(), 23u);

    // Verify some values (ground truth from oracle)
    EXPECT_NEAR(atom_radii[0].second, 0.77, 1e-4);
    EXPECT_NEAR(atom_radii[1].second, 0.65, 1e-4);
    EXPECT_NEAR(atom_radii[2].second, 0.6, 1e-4);
    EXPECT_NEAR(atom_radii[3].second, 0.35, 1e-4);
    EXPECT_NEAR(atom_radii[4].second, 1.1, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}