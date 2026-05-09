#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetResidueAveragePosition) {
    std::string pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
    gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

    // Test 1: Valid CID for residue A/10
    std::string cid1 = "//A/10";
    std::vector<double> avg1 = coot::molecule_t::get_residue_average_position_gemmi(st, cid1);
    EXPECT_EQ(avg1.size(), 3u);
    if (avg1.size() == 3) {
        EXPECT_NEAR(avg1[0], 27.2126, 1e-3);
        EXPECT_NEAR(avg1[1], 11.192, 1e-3);
        EXPECT_NEAR(avg1[2], 67.5274, 1e-3);
    }

    // Test 2: Valid CID for residue A/50
    std::string cid2 = "//A/50";
    std::vector<double> avg2 = coot::molecule_t::get_residue_average_position_gemmi(st, cid2);
    EXPECT_EQ(avg2.size(), 3u);
    if (avg2.size() == 3) {
        EXPECT_NEAR(avg2[0], 6.15482, 1e-3);
        EXPECT_NEAR(avg2[1], -3.97091, 1e-3);
        EXPECT_NEAR(avg2[2], 56.4205, 1e-3);
    }

    // Test 3: Invalid CID (nonexistent residue)
    std::string cid3 = "//B/999";
    std::vector<double> avg3 = coot::molecule_t::get_residue_average_position_gemmi(st, cid3);
    EXPECT_EQ(avg3.size(), 0u);

    // Test 4: Valid CID for first residue A/1
    std::string cid4 = "//A/1";
    std::vector<double> avg4 = coot::molecule_t::get_residue_average_position_gemmi(st, cid4);
    EXPECT_EQ(avg4.size(), 3u);
    if (avg4.size() == 3) {
        EXPECT_NEAR(avg4[0], 13.8456, 1e-3);
        EXPECT_NEAR(avg4[1], 5.18788, 1e-3);
        EXPECT_NEAR(avg4[2], 75.0207, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}