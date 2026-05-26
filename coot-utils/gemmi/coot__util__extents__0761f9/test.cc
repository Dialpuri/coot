#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ExtentsFunction) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test case 1: All atoms (empty chain_id means no filter)
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result1 = coot::util::extents_gemmi(st);
    
    EXPECT_NEAR(result1.first.x(), 1.593, 1e-3);
    EXPECT_NEAR(result1.first.y(), -14.67, 1e-3);
    EXPECT_NEAR(result1.first.z(), 20.067, 1e-3);
    EXPECT_NEAR(result1.second.x(), 51.369, 1e-3);
    EXPECT_NEAR(result1.second.y(), 25.796, 1e-3);
    EXPECT_NEAR(result1.second.z(), 82.943, 1e-3);
    
    // Test case 2: Specific chain A
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result2 = coot::util::extents_gemmi(st, "A");
    
    EXPECT_NEAR(result2.first.x(), 1.593, 1e-3);
    EXPECT_NEAR(result2.first.y(), -14.67, 1e-3);
    EXPECT_NEAR(result2.first.z(), 20.067, 1e-3);
    EXPECT_NEAR(result2.second.x(), 51.369, 1e-3);
    EXPECT_NEAR(result2.second.y(), 25.796, 1e-3);
    EXPECT_NEAR(result2.second.z(), 82.943, 1e-3);
    
    // Test case 3: Select a specific residue (using "*" wildcard for all chains, equivalent to all atoms)
    std::pair<clipper::Coord_orth, clipper::Coord_orth> result3 = coot::util::extents_gemmi(st, "*");
    
    EXPECT_NEAR(result3.first.x(), 1.593, 1e-3);
    EXPECT_NEAR(result3.first.y(), -14.67, 1e-3);
    EXPECT_NEAR(result3.first.z(), 20.067, 1e-3);
    EXPECT_NEAR(result3.second.x(), 51.369, 1e-3);
    EXPECT_NEAR(result3.second.y(), 25.796, 1e-3);
    EXPECT_NEAR(result3.second.z(), 82.943, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}