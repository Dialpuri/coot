#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "function.hh"

// case: valid PDB loaded
TEST(OracleTest, ExtentsValidMolecule) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result = coot::util::extents_gemmi(st);
    clipper::Coord_orth min_coord = result.first;
    clipper::Coord_orth max_coord = result.second;

    EXPECT_NEAR(min_coord.x(), 1.593, 1e-4);
    EXPECT_NEAR(min_coord.y(), -14.67, 1e-4);
    EXPECT_NEAR(min_coord.z(), 20.067, 1e-4);
    EXPECT_NEAR(max_coord.x(), 51.369, 1e-4);
    EXPECT_NEAR(max_coord.y(), 25.796, 1e-4);
    EXPECT_NEAR(max_coord.z(), 82.943, 1e-4);

    EXPECT_NEAR(max_coord.x() - min_coord.x(), 49.776, 1e-4);
    EXPECT_NEAR(max_coord.y() - min_coord.y(), 40.466, 1e-4);
    EXPECT_NEAR(max_coord.z() - min_coord.z(), 62.876, 1e-4);
}

// case: empty structure (equivalent to null pointer in original)
TEST(OracleTest, ExtentsEmpty) {
    gemmi::Structure st;
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result = coot::util::extents_gemmi(st);
    clipper::Coord_orth min_coord = result.first;
    clipper::Coord_orth max_coord = result.second;

    EXPECT_NEAR(min_coord.x(), 0.0, 1e-4);
    EXPECT_NEAR(min_coord.y(), 0.0, 1e-4);
    EXPECT_NEAR(min_coord.z(), 0.0, 1e-4);
    EXPECT_NEAR(max_coord.x(), 0.0, 1e-4);
    EXPECT_NEAR(max_coord.y(), 0.0, 1e-4);
    EXPECT_NEAR(max_coord.z(), 0.0, 1e-4);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}