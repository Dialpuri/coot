#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "function.hh"

TEST(OracleTest, shift_gemmi) {
    // Load a PDB file into a gemmi::Structure
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get first atom for before/after comparison
    gemmi::Model &model = st.models[0];
    gemmi::Chain &chain = model.chains[0];
    gemmi::Residue &res = chain.residues[0];
    gemmi::Atom &at = res.atoms[0];

    // Verify atom name matches gemmi (unpadded) representation
    EXPECT_EQ(at.name, "N");

    // Case 1: non-zero shift — verifies the mutation path
    {
        double x_before = at.pos.x;
        double y_before = at.pos.y;
        double z_before = at.pos.z;
        EXPECT_NEAR(x_before, 12.334, 1e-4);
        EXPECT_NEAR(y_before, 2.772, 1e-4);
        EXPECT_NEAR(z_before, 74.951, 1e-4);

        clipper::Coord_orth shift(0.5, 1.0, -0.25);
        coot::util::shift_gemmi(st, shift);

        EXPECT_NEAR(at.pos.x, 12.834, 1e-4);
        EXPECT_NEAR(at.pos.y, 3.772, 1e-4);
        EXPECT_NEAR(at.pos.z, 74.701, 1e-4);

        // Verify coordinates were actually modified
        EXPECT_FALSE(at.pos.x == x_before && at.pos.y == y_before && at.pos.z == z_before);
    }

    // Case 2: zero shift — verifies guard path (coordinates unchanged)
    {
        double x_before = at.pos.x;
        double y_before = at.pos.y;
        double z_before = at.pos.z;

        clipper::Coord_orth shift(0.0, 0.0, 0.0);
        coot::util::shift_gemmi(st, shift);

        EXPECT_NEAR(at.pos.x, 12.834, 1e-4);
        EXPECT_NEAR(at.pos.y, 3.772, 1e-4);
        EXPECT_NEAR(at.pos.z, 74.701, 1e-4);

        // Verify coordinates are unchanged (zero shift is a no-op)
        EXPECT_DOUBLE_EQ(at.pos.x, x_before);
        EXPECT_DOUBLE_EQ(at.pos.y, y_before);
        EXPECT_DOUBLE_EQ(at.pos.z, z_before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}