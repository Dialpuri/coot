#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include "function.hh"

TEST(OracleTest, atom_constructor) {
    // Case 1: typical atom construction via the gemmi::Position constructor
    {
        std::string atom_name = "CA";
        std::string ele = "C";
        gemmi::Position pos(10.5f, 20.3f, 30.7f);
        std::string altloc = "";
        float dbf = 25.0f;

        coot::minimol::atom at = coot::minimol::atom_gemmi(atom_name, ele, pos, altloc, dbf);

        EXPECT_EQ(at.name, "CA");
        EXPECT_EQ(at.element, "C");
        EXPECT_NEAR(at.pos.x, 10.5f, 1e-4);
        EXPECT_NEAR(at.pos.y, 20.3f, 1e-4);
        EXPECT_NEAR(at.pos.z, 30.7f, 1e-4);
        EXPECT_TRUE(at.altLoc.empty());
        EXPECT_NEAR(at.occupancy, 1.0f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 25.0f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }

    // Case 2: atom with alternate location and different dbf
    {
        std::string atom_name = "N";
        std::string ele = "N";
        gemmi::Position pos(1.0f, 2.0f, 3.0f);
        std::string altloc = "A";
        float dbf = 15.0f;

        coot::minimol::atom at = coot::minimol::atom_gemmi(atom_name, ele, pos, altloc, dbf);

        EXPECT_EQ(at.name, "N");
        EXPECT_EQ(at.element, "N");
        EXPECT_NEAR(at.pos.x, 1.0f, 1e-4);
        EXPECT_NEAR(at.pos.y, 2.0f, 1e-4);
        EXPECT_NEAR(at.pos.z, 3.0f, 1e-4);
        EXPECT_EQ(at.altLoc, "A");
        EXPECT_NEAR(at.occupancy, 1.0f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 15.0f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }

    // Case 3: edge case - zero position, zero dbf
    {
        std::string atom_name = "OXT";
        std::string ele = "O";
        gemmi::Position pos(0.0f, 0.0f, 0.0f);
        std::string altloc = "";
        float dbf = 0.0f;

        coot::minimol::atom at = coot::minimol::atom_gemmi(atom_name, ele, pos, altloc, dbf);

        EXPECT_EQ(at.name, "OXT");
        EXPECT_EQ(at.element, "O");
        EXPECT_NEAR(at.pos.x, 0.0f, 1e-4);
        EXPECT_NEAR(at.pos.y, 0.0f, 1e-4);
        EXPECT_NEAR(at.pos.z, 0.0f, 1e-4);
        EXPECT_TRUE(at.altLoc.empty());
        EXPECT_NEAR(at.occupancy, 1.0f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 0.0f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
