#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, read_file) {
    // Case 1: valid PDB
    {
        coot::minimol::molecule m;
        std::string pdb = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        int ierr = coot::minimol::read_file_gemmi(m, pdb);
        EXPECT_EQ(ierr, 0);
        EXPECT_EQ(m.fragments.size(), 2);
        EXPECT_TRUE(m.has_atoms());
    }

    // Case 2: non-existent file
    {
        coot::minimol::molecule m2;
        std::string bad_pdb = "/nonexistent/path/file.pdb";
        int ierr2 = coot::minimol::read_file_gemmi(m2, bad_pdb);
        EXPECT_EQ(ierr2, 15);
        EXPECT_EQ(m2.fragments.size(), 0);
        EXPECT_FALSE(m2.has_atoms());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}