#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, GetFirstResidueInChain) {
    // case 1: valid chain — should return first residue
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        ASSERT_FALSE(st.models.empty()) << "no model";

        const gemmi::Chain& chain = st.models[0].chains[0];
        const gemmi::Residue* first_res = coot::util::get_first_residue_in_chain_gemmi(chain);

        EXPECT_TRUE(first_res != nullptr) << "first_residue should not be nullptr for chain " << chain.name;
        if (first_res) {
            EXPECT_EQ(first_res->seqid.num.value, 1);
            EXPECT_EQ(first_res->name, "MET");
        }
    }

    // case 2: empty chain — should return nullptr
    {
        gemmi::Chain empty_chain;
        empty_chain.name = "A";

        const gemmi::Residue* first_res = coot::util::get_first_residue_in_chain_gemmi(empty_chain);

        EXPECT_FALSE(first_res != nullptr) << "first_residue should be nullptr for empty chain";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}