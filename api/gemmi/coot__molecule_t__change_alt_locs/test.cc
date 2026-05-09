#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ChangeAltLocs) {
    // Load the test molecule with explicit hydrogens
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Test with chain A, residue 1 using different modes
    std::string cid = "//A/1";

    // Test 1: residue mode
    {
        std::string change_mode = "residue";
        int status = coot::molecule_t::change_alt_locs_gemmi(st, cid, change_mode);
        EXPECT_EQ(status, 0);
    }

    // Test 2: side-chain mode
    {
        std::string change_mode2 = "side-chain";
        int status = coot::molecule_t::change_alt_locs_gemmi(st, cid, change_mode2);
        EXPECT_EQ(status, 0);
    }

    // Test 3: main-chain mode
    {
        std::string change_mode3 = "main-chain";
        int status = coot::molecule_t::change_alt_locs_gemmi(st, cid, change_mode3);
        EXPECT_EQ(status, 0);
    }

    // Test 4: atom list mode (CB)
    {
        std::string change_mode4 = "CB";
        int status = coot::molecule_t::change_alt_locs_gemmi(st, cid, change_mode4);
        EXPECT_EQ(status, 0);
    }

    // Test 5: residue mode on residue 10
    {
        std::string cid2 = "//A/10";
        std::string change_mode = "residue";
        int status = coot::molecule_t::change_alt_locs_gemmi(st, cid2, change_mode);
        EXPECT_EQ(status, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}