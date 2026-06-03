#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

TEST(OracleTest, closest_approach_gemmi) {
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    const gemmi::SpaceGroup *sg = st.find_spacegroup();
    int nsymm = (sg) ? static_cast<int>(sg->operations().sym_ops.size()) : 0;
    EXPECT_EQ(nsymm, 4);

    // Case 1: valid points
    {
        clipper::Coord_orth moving_point(10.0, 20.0, 30.0);
        clipper::Coord_orth reference_point(0.0, 0.0, 0.0);

        clipper::Coord_orth result = closest_approach_gemmi(moving_point, reference_point, st);

        EXPECT_NEAR(result.x(), -10.0, 1e-4);
        EXPECT_NEAR(result.y(), -15.6555, 1e-4);
        EXPECT_NEAR(result.z(), 6.0505, 1e-4);
    }

    // Case 2: same point for moving and reference
    {
        clipper::Coord_orth moving_point(10.0, 20.0, 30.0);
        clipper::Coord_orth reference_point(10.0, 20.0, 30.0);

        clipper::Coord_orth result = closest_approach_gemmi(moving_point, reference_point, st);

        EXPECT_NEAR(result.x(), 10.0, 1e-4);
        EXPECT_NEAR(result.y(), 20.0, 1e-4);
        EXPECT_NEAR(result.z(), 30.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}