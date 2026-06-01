#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

TEST(OracleTest, gemmi_manager) {
    // Case: basic structure construction (empty)
    {
        gemmi::Structure st;
        EXPECT_EQ(st.models.size(), 0);
        EXPECT_EQ(gemmi::count_atom_sites(st), 0);
    }

    // Case: load PDB and verify content
    {
        gemmi::Structure st = CXXCreator::CXXCreator_gemmi("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st.models.size(), 1);
        EXPECT_GT(gemmi::count_atom_sites(st), 0);

        // Get the first model (0-based indexing in gemmi)
        EXPECT_FALSE(st.models.empty());
        EXPECT_GT(st.models[0].chains.size(), 0u);
    }

    // Case: loading same PDB twice produces independent structures (gemmi reads are independent calls)
    {
        gemmi::Structure st1 = CXXCreator::CXXCreator_gemmi("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Structure st2 = CXXCreator::CXXCreator_gemmi("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        EXPECT_EQ(st1.models.size(), 1);
        EXPECT_GT(gemmi::count_atom_sites(st1), 0);
        EXPECT_EQ(st2.models.size(), 1);
        EXPECT_GT(gemmi::count_atom_sites(st2), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}