#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, chains_in_atom_selection) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Case 1: select chain A
    {
        std::string atom_selection = "//A";
        std::vector<std::string> chains = coot::util::chains_in_atom_selection_gemmi(model, atom_selection);
        EXPECT_EQ(chains.size(), 1u);
        EXPECT_EQ(chains[0], "A");
    }

    // Case 2: select chain B (water)
    {
        std::string atom_selection = "//B";
        std::vector<std::string> chains = coot::util::chains_in_atom_selection_gemmi(model, atom_selection);
        EXPECT_EQ(chains.size(), 1u);
        EXPECT_EQ(chains[0], "B");
    }

    // Case 3: invalid chain selection
    {
        std::string atom_selection = "//Z";
        std::vector<std::string> chains = coot::util::chains_in_atom_selection_gemmi(model, atom_selection);
        EXPECT_EQ(chains.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}