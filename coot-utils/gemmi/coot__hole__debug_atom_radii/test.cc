#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, coot_hole_debug_atom_radii_gemmi) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1);

        coot::protein_geometry geometry;
        geometry.init_standard();

        auto radii = coot::hole::assign_vdw_radii_gemmi(st.models[0], geometry);

        // debug_atom_radii_gemmi is void and prints to stdout — verify it doesn't throw
        EXPECT_NO_THROW(coot::hole::debug_atom_radii_gemmi(st.models[0], radii));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}