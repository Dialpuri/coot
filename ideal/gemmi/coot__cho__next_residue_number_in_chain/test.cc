#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, next_residue_number_in_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain_a = model.chains[0];

    // Case 1: valid chain, new_res_no_by_hundreds = false
    {
        bool by_hundreds = false;
        std::pair<short int, int> result = coot::cho::next_residue_number_in_chain_gemmi(&chain_a, by_hundreds);
        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, 299);
    }

    // Case 2: valid chain, new_res_no_by_hundreds = true
    {
        bool by_hundreds = true;
        std::pair<short int, int> result = coot::cho::next_residue_number_in_chain_gemmi(&chain_a, by_hundreds);
        EXPECT_EQ(result.first, 1);
        EXPECT_EQ(result.second, 301);
    }

    // Case 3: null chain — edge case
    {
        bool by_hundreds = false;
        const gemmi::Chain* null_chain = nullptr;
        std::pair<short int, int> result = coot::cho::next_residue_number_in_chain_gemmi(null_chain, by_hundreds);
        EXPECT_EQ(result.first, 0);
        EXPECT_EQ(result.second, 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}