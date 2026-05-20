#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_make_atom) {
    // case 1: typical atom with all fields set
    {
        coot::minimol::atom a;
        a.name = "CA";
        a.element = "C";
        a.pos = {1.5f, 2.3f, 3.7f};
        a.occupancy = 1.0f;
        a.temperature_factor = 20.5f;

        gemmi::Atom at = coot::minimol::make_atom_gemmi(a);

        EXPECT_STREQ(at.name.c_str(), "CA");
        EXPECT_STREQ(at.element.name(), "C");
        EXPECT_NEAR(at.pos.x, 1.5, 1e-4);
        EXPECT_NEAR(at.pos.y, 2.3, 1e-4);
        EXPECT_NEAR(at.pos.z, 3.7, 1e-4);
        EXPECT_NEAR(at.occ, 1.0, 1e-4);
        EXPECT_NEAR(at.b_iso, 20.5, 1e-4);
    }

    // case 2: minimal atom (defaults for pos/occupancy/temp)
    {
        coot::minimol::atom a;
        a.name = "N";
        a.element = "N";

        gemmi::Atom at = coot::minimol::make_atom_gemmi(a);

        EXPECT_STREQ(at.name.c_str(), "N");
        EXPECT_STREQ(at.element.name(), "N");
        EXPECT_NEAR(at.pos.x, 0.0, 1e-4);
        EXPECT_NEAR(at.pos.y, 0.0, 1e-4);
        EXPECT_NEAR(at.pos.z, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}