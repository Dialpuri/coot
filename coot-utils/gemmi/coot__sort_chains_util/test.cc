#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, sort_chains_util) {
    // Setup: load PDB to get real chain data
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: sort chain IDs using sort_chains_util_gemmi as comparator
    {
        std::vector<std::string> chain_ids;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                chain_ids.push_back(chain.name);
            }
        }

        std::sort(chain_ids.begin(), chain_ids.end(), coot::sort_chains_util_gemmi);

        EXPECT_EQ(chain_ids.size(), 2);
        EXPECT_EQ(chain_ids[0], "A");
        EXPECT_EQ(chain_ids[1], "B");
    }

    // Case 2: direct comparison A vs B
    {
        bool r = coot::sort_chains_util_gemmi("A", "B");
        EXPECT_TRUE(r);
    }

    // Case 3: direct comparison B vs A (reverse order)
    {
        bool r = coot::sort_chains_util_gemmi("B", "A");
        EXPECT_FALSE(r);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}