#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_last_residue_in_chain) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "Model 1 should exist";

    // case 1: chain "A" exists and has residues
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Chain& chain = model.chains[0]; // chain A
        const gemmi::Residue* last_res = coot::util::get_last_residue_in_chain_gemmi(chain);
        EXPECT_TRUE(last_res != nullptr) << "last_residue should be found for valid chain";
        if (last_res) {
            EXPECT_EQ(last_res->seqid.num.value, 298);
            EXPECT_EQ(last_res->name, "LEU");
        }
    }

    // case 2: non-existent chain (guard clause)
    {
        const gemmi::Model& model = st.models[0];
        if (model.chains.size() > 99) {
            const gemmi::Chain& chain = model.chains[99];
            const gemmi::Residue* last_res = coot::util::get_last_residue_in_chain_gemmi(chain);
            EXPECT_TRUE(last_res != nullptr);
        } else {
            // Simulate out-of-range by using a dummy chain with no residues
            gemmi::Chain dummy_chain;
            dummy_chain.name = "X";
            const gemmi::Residue* last_res = coot::util::get_last_residue_in_chain_gemmi(dummy_chain);
            EXPECT_EQ(last_res, nullptr) << "should return nullptr for empty chain";
        }
    }

    // case 3: null chain pointer — cannot test with references; skip
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}