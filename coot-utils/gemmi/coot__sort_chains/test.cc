#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, sort_chains) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // case: number of models
    int nModels = (int)st.models.size();
    EXPECT_EQ(nModels, 1);

    // case: chains before sort_chains
    gemmi::Model &model = st.models[0];
    int nChains = (int)model.chains.size();
    EXPECT_EQ(nChains, 2);
    EXPECT_STREQ(model.chains[0].name.c_str(), "A");
    EXPECT_STREQ(model.chains[1].name.c_str(), "B");

    // case: call sort_chains_gemmi
    coot::sort_chains_gemmi(st);

    // case: chains after sort_chains
    nChains = (int)model.chains.size();
    EXPECT_EQ(nChains, 2);
    EXPECT_STREQ(model.chains[0].name.c_str(), "A");
    EXPECT_STREQ(model.chains[1].name.c_str(), "B");
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}