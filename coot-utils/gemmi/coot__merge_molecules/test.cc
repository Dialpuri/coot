#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, merge_molecules) {
    // case 1: merge two distinct molecules (same PDB)
    {
        gemmi::Structure st1 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Structure st2 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int chains_before = st1.models[0].chains.size();

        std::vector<gemmi::Structure> others;
        others.push_back(st2);

        int mol2_chains = st2.models[0].chains.size();

        coot::merge_molecules_gemmi(st1, others);

        int chains_after = st1.models[0].chains.size();
        EXPECT_EQ(chains_before, 2);
        EXPECT_EQ(mol2_chains, 2);
        EXPECT_EQ(chains_after, 4);
    }

    // case 2: merge with empty others vector
    {
        gemmi::Structure st1 = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int chains_before = st1.models[0].chains.size();

        std::vector<gemmi::Structure> others;

        coot::merge_molecules_gemmi(st1, others);

        int chains_after = st1.models[0].chains.size();
        EXPECT_EQ(chains_before, 2);
        EXPECT_EQ(chains_after, 2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}