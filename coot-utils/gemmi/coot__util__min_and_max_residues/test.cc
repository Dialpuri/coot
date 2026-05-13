#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, min_and_max_residues) {
    // Case: valid chain with residues
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_GT(st.models.size(), 0);
        gemmi::Model& model = st.models[0];
        ASSERT_GT(model.chains.size(), 0);
        gemmi::Chain& chain = model.chains[0];

        std::pair<int, int> result = coot::util::min_and_max_residues_gemmi(chain);
        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, 298);
    }

    // Case: null chain pointer (edge case)
    {
        std::pair<int, int> result = coot::util::min_and_max_residues_gemmi(static_cast<gemmi::Chain*>(nullptr));
        EXPECT_EQ(result.first, 9999);
        EXPECT_EQ(result.second, -9999);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}