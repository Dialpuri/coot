#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, NumberOfChains) {
    // case 1: valid structure with 2 chains
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        int nchains = coot::util::number_of_chains_gemmi(&st);
        EXPECT_EQ(nchains, 2);
    }

    // case 2: null pointer (badness case)
    {
        gemmi::Structure* null_st = nullptr;
        int nchains = coot::util::number_of_chains_gemmi(null_st);
        EXPECT_EQ(nchains, -1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}