#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_ncs_related_chains_example_pdb) {
    // case: setup
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1u) << "PDB load failed";

    // case: standard structure with 2 chains A and B
    {
        const gemmi::Model& model = st.models[0];
        std::vector<std::vector<const gemmi::Chain*>> ncs_chains = coot::ncs_related_chains_gemmi(model);
        EXPECT_EQ(ncs_chains.size(), 2u);
        for (const auto& group : ncs_chains) {
            EXPECT_EQ(group.size(), 1u);
            for (const auto& chain : group) {
                EXPECT_TRUE(chain->name == "A" || chain->name == "B");
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}