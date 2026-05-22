#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_second_residue_in_chain) {
    // Setup: load the example PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    const gemmi::Model* model = &st.models[0];
    
    // case: chain with multiple residues (should return second residue)
    {
        const gemmi::Chain& chain = model->chains[0];
        const gemmi::Residue* second_res = coot::util::get_second_residue_in_chain_gemmi(chain);
        EXPECT_NE(second_res, nullptr);
        if (second_res) {
            EXPECT_EQ(second_res->seqid.num.value, 2);
            EXPECT_EQ(second_res->name, "GLU");
        }
    }

    // case: null chain pointer (should return nullptr)
    {
        const gemmi::Residue* second_res = coot::util::get_second_residue_in_chain_gemmi(nullptr, "A");
        EXPECT_EQ(second_res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}