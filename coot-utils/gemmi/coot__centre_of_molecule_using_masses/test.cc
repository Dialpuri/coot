#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, centre_of_molecule_using_masses) {
    // Case 1: Valid molecule with atoms
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        auto result = coot::centre_of_molecule_using_masses_gemmi(st);
        EXPECT_TRUE(result.first);
        if (result.first) {
            clipper::Coord_orth coord = result.second;
            EXPECT_NEAR(coord.x(), 29.438, 1e-4);
            EXPECT_NEAR(coord.y(), 5.05099, 1e-4);
            EXPECT_NEAR(coord.z(), 50.2175, 1e-4);
        }
    }

    // Case 2: Empty molecule (no atoms)
    {
        gemmi::Structure empty_st;
        auto result = coot::centre_of_molecule_using_masses_gemmi(empty_st);
        EXPECT_FALSE(result.first);
        if (result.first) {
            clipper::Coord_orth coord = result.second;
            EXPECT_NEAR(coord.x(), 0, 1e-4);
            EXPECT_NEAR(coord.y(), 0, 1e-4);
            EXPECT_NEAR(coord.z(), 0, 1e-4);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}