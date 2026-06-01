#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, residues_in_order_p_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Case 1: chain A (ordered protein residues)
    {
        const gemmi::Chain* chain_p = model.find_chain("A");
        EXPECT_NE(chain_p, nullptr);
        EXPECT_EQ(static_cast<int>(chain_p->residues.size()), 267);
        bool result = coot::residues_in_order_p_gemmi(chain_p);
        EXPECT_EQ(result, true);
    }

    // Case 2: null chain pointer (guarded path — returns true)
    {
        const gemmi::Chain* chain_p = nullptr;
        bool result = coot::residues_in_order_p_gemmi(chain_p);
        EXPECT_EQ(result, true);
    }

    // Case 3: chain B (1 HOH residue — trivially ordered)
    {
        const gemmi::Chain* chain_p = model.find_chain("B");
        EXPECT_NE(chain_p, nullptr);
        EXPECT_EQ(static_cast<int>(chain_p->residues.size()), 1);
        bool result = coot::residues_in_order_p_gemmi(chain_p);
        EXPECT_EQ(result, true);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}