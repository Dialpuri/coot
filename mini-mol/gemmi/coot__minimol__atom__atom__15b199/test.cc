#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, minimol_atom) {

    // Case 1: typical atom with standard values
    {
        std::string atom_name = "CA";
        std::string ele = "C";
        float x = 1.23f;
        float y = 4.56f;
        float z = 7.89f;
        std::string altloc = "";
        float occupancy = 1.0f;
        float dbf = 20.0f;

        coot::minimol::atom at = coot::minimol::atom::atom_gemmi(atom_name, ele, x, y, z, altloc, occupancy, dbf);

        EXPECT_EQ(at.name, "CA");
        EXPECT_EQ(at.element, "C");
        EXPECT_NEAR(at.pos.x, 1.23f, 1e-4);
        EXPECT_NEAR(at.pos.y, 4.56f, 1e-4);
        EXPECT_NEAR(at.pos.z, 7.89f, 1e-4);
        EXPECT_EQ(at.altLoc, "");
        EXPECT_NEAR(at.occupancy, 1.0f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 20.0f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }

    // Case 2: atom with alternate location and partial occupancy
    {
        std::string atom_name = "N";
        std::string ele = "N";
        float x = -3.14f;
        float y = 0.0f;
        float z = 15.5f;
        std::string altloc = "A";
        float occupancy = 0.5f;
        float dbf = 35.7f;

        coot::minimol::atom at = coot::minimol::atom::atom_gemmi(atom_name, ele, x, y, z, altloc, occupancy, dbf);

        EXPECT_EQ(at.name, "N");
        EXPECT_EQ(at.element, "N");
        EXPECT_NEAR(at.pos.x, -3.14f, 1e-4);
        EXPECT_NEAR(at.pos.y, 0.0f, 1e-4);
        EXPECT_NEAR(at.pos.z, 15.5f, 1e-4);
        EXPECT_EQ(at.altLoc, "A");
        EXPECT_NEAR(at.occupancy, 0.5f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 35.7f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }

    // Case 3: atom with zero coordinates and no occupancy/B-factor
    {
        std::string atom_name = "OXT";
        std::string ele = "O";
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
        std::string altloc = "";
        float occupancy = 0.0f;
        float dbf = 0.0f;

        coot::minimol::atom at = coot::minimol::atom::atom_gemmi(atom_name, ele, x, y, z, altloc, occupancy, dbf);

        EXPECT_EQ(at.name, "OXT");
        EXPECT_EQ(at.element, "O");
        EXPECT_NEAR(at.pos.x, 0.0f, 1e-4);
        EXPECT_NEAR(at.pos.y, 0.0f, 1e-4);
        EXPECT_NEAR(at.pos.z, 0.0f, 1e-4);
        EXPECT_EQ(at.altLoc, "");
        EXPECT_NEAR(at.occupancy, 0.0f, 1e-4);
        EXPECT_NEAR(at.temperature_factor, 0.0f, 1e-4);
        EXPECT_EQ(at.int_user_data, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
