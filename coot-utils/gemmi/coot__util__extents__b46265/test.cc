#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include "function.hh"

TEST(OracleTest, Extents) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    std::pair<clipper::Coord_orth, clipper::Coord_orth> e = coot::util::extents_gemmi(st);
    EXPECT_NEAR(e.first.x(), 1.593, 1e-3);
    EXPECT_NEAR(e.first.y(), -14.67, 1e-3);
    EXPECT_NEAR(e.first.z(), 20.067, 1e-3);
    EXPECT_NEAR(e.second.x(), 51.369, 1e-3);
    EXPECT_NEAR(e.second.y(), 25.796, 1e-3);
    EXPECT_NEAR(e.second.z(), 82.943, 1e-3);
    
    clipper::Coord_orth range = e.second - e.first;
    EXPECT_NEAR(range.x(), 49.776, 1e-3);
    EXPECT_NEAR(range.y(), 40.466, 1e-3);
    EXPECT_NEAR(range.z(), 62.876, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}