#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

gemmi::Structure read_pdb(const std::string &filename) {
    return gemmi::read_pdb_file(filename);
}

TEST(OracleTest, get_residue_mid_point) {
    // case 1: valid residue spec - residue 10 in chain A
    {
        gemmi::Structure st = read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        std::pair<bool, clipper::Coord_orth> result = 
            coot::util::get_residue_mid_point_gemmi("A", 10, "", st);
        EXPECT_TRUE(result.first);
        if (result.first) {
            EXPECT_NEAR(result.second.x(), 27.2126, 1e-4);
            EXPECT_NEAR(result.second.y(), 11.192, 1e-4);
            EXPECT_NEAR(result.second.z(), 67.5274, 1e-4);
        }
    }

    // case 2: invalid residue spec - non-existent residue
    {
        gemmi::Structure st = read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        std::pair<bool, clipper::Coord_orth> result = 
            coot::util::get_residue_mid_point_gemmi("A", 9999, "", st);
        EXPECT_FALSE(result.first);
    }

    // case 3: empty structure
    {
        gemmi::Structure st;
        std::pair<bool, clipper::Coord_orth> result = 
            coot::util::get_residue_mid_point_gemmi("A", 10, "", st);
        EXPECT_FALSE(result.first);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}