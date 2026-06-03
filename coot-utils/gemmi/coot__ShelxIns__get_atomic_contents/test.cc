#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_atomic_contents_gemmi) {
    // case: main
    {
        // Load PDB into gemmi::Structure
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        // Create ShelxIns and call get_atomic_contents_gemmi
        coot::ShelxIns shelx;
        std::map<std::string, unsigned int> contents = shelx.get_atomic_contents_gemmi(st);

        EXPECT_EQ(contents.size(), 4u);
        EXPECT_EQ(contents[" C"], 1364u);
        EXPECT_EQ(contents[" N"], 360u);
        EXPECT_EQ(contents[" O"], 375u);
        EXPECT_EQ(contents[" S"], 8u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}