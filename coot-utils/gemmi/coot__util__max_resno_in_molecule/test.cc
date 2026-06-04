#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, max_resno_in_molecule) {
    // Case 1: standard protein PDB
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1);

        std::pair<bool, int> result = coot::util::max_resno_in_molecule_gemmi(st);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second, 298);
    }

    // Case 2: protein with hydrogens
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        EXPECT_EQ(st.models.size(), 1);

        std::pair<bool, int> result = coot::util::max_resno_in_molecule_gemmi(st);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second, 298);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}