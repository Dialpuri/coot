#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, minimol_atom) {
  // Case 1: typical atom (CA from residue 10, chain A)
  {
    std::string atom_name = "CA";
    std::string element   = "C";
    gemmi::Vec3 pos(10.5f, 20.3f, 30.7f);
    std::string altloc    = "";
    float occupancy       = 1.0f;
    float dbf             = 15.5f;

    coot::minimol::atom a = coot::minimol::atom::atom_gemmi(atom_name, element, pos, altloc, occupancy, dbf);

    EXPECT_EQ(a.name, "CA");
    EXPECT_EQ(a.element, "C");
    EXPECT_NEAR(a.pos.x, 10.5f, 1e-4);
    EXPECT_NEAR(a.pos.y, 20.3f, 1e-4);
    EXPECT_NEAR(a.pos.z, 30.7f, 1e-4);
    EXPECT_EQ(a.altLoc, "");
    EXPECT_NEAR(a.occupancy, 1.0f, 1e-4);
    EXPECT_NEAR(a.temperature_factor, 15.5f, 1e-4);
    EXPECT_EQ(a.int_user_data, -1);
  }

  // Case 2: atom with alternate location and partial occupancy
  {
    std::string atom_name = "OG";
    std::string element   = "O";
    gemmi::Vec3 pos(5.1f, 5.2f, 5.3f);
    std::string altloc    = "A";
    float occupancy       = 0.6f;
    float dbf             = 22.0f;

    coot::minimol::atom a = coot::minimol::atom::atom_gemmi(atom_name, element, pos, altloc, occupancy, dbf);

    EXPECT_EQ(a.name, "OG");
    EXPECT_EQ(a.element, "O");
    EXPECT_NEAR(a.pos.x, 5.1f, 1e-4);
    EXPECT_NEAR(a.pos.y, 5.2f, 1e-4);
    EXPECT_NEAR(a.pos.z, 5.3f, 1e-4);
    EXPECT_EQ(a.altLoc, "A");
    EXPECT_NEAR(a.occupancy, 0.6f, 1e-4);
    EXPECT_NEAR(a.temperature_factor, 22.0f, 1e-4);
    EXPECT_EQ(a.int_user_data, -1);
  }

  // Case 3: edge case — zero occupancy, zero B-factor, negative coordinates
  {
    std::string atom_name = "H1";
    std::string element   = "H";
    gemmi::Vec3 pos(-1.0f, 0.0f, 0.0f);
    std::string altloc    = "";
    float occupancy       = 0.0f;
    float dbf             = 0.0f;

    coot::minimol::atom a = coot::minimol::atom::atom_gemmi(atom_name, element, pos, altloc, occupancy, dbf);

    EXPECT_EQ(a.name, "H1");
    EXPECT_EQ(a.element, "H");
    EXPECT_NEAR(a.pos.x, -1.0f, 1e-4);
    EXPECT_NEAR(a.pos.y, 0.0f, 1e-4);
    EXPECT_NEAR(a.pos.z, 0.0f, 1e-4);
    EXPECT_EQ(a.altLoc, "");
    EXPECT_NEAR(a.occupancy, 0.0f, 1e-4);
    EXPECT_NEAR(a.temperature_factor, 0.0f, 1e-4);
    EXPECT_EQ(a.int_user_data, -1);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
