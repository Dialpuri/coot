#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_default_constructor) {
    // case: default
    coot::minimol::atom default_atom;

    // int_user_data
    EXPECT_EQ(default_atom.int_user_data, -1);
    // occupancy
    EXPECT_NEAR(default_atom.occupancy, -1, 1e-4);
    // temperature_factor
    EXPECT_NEAR(default_atom.temperature_factor, -1, 1e-4);
    // name (empty string)
    EXPECT_EQ(default_atom.name, "");
    // element (empty string)
    EXPECT_EQ(default_atom.element, "");
    // altLoc (empty string)
    EXPECT_EQ(default_atom.altLoc, "");
    // pos (default constructed clipper::Coord_orth - near zero)
    EXPECT_NEAR(default_atom.pos.x, 0.0, 1e-4);
    EXPECT_NEAR(default_atom.pos.y, 0.0, 1e-4);
    EXPECT_NEAR(default_atom.pos.z, 0.0, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}