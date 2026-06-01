#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static const gemmi::Residue* find_residue(
    const gemmi::Model& model, const std::string& chain_id, int seq_num)
{
    for (const auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (const auto& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, peptide_C_and_N_are_close_p) {
    gemmi::Structure st =
        gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    const gemmi::Model& model = st.models[0];

    // Get two adjacent residues (A/1 MET and A/2 GLU)
    const gemmi::Residue* res1 = find_residue(model, "A", 1);
    const gemmi::Residue* res2 = find_residue(model, "A", 2);
    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    // Get a distant residue (A/100 PRO)
    const gemmi::Residue* res_distant = find_residue(model, "A", 100);

    // Case 1: Adjacent residues A/1 (MET) and A/2 (GLU) — should be close
    {
        std::pair<bool, bool> result =
            coot::peptide_C_and_N_are_close_p_gemmi(res1, res2);
        EXPECT_TRUE(result.first);   // are_close
        EXPECT_FALSE(result.second); // order_switch
    }

    // Case 2: Distant residues A/1 (MET) and A/100 (PRO) — should NOT be close
    {
        if (res_distant) {
            std::pair<bool, bool> result =
                coot::peptide_C_and_N_are_close_p_gemmi(res1, res_distant);
            EXPECT_FALSE(result.first);  // are_close
            EXPECT_FALSE(result.second); // order_switch
        }
    }

    // Case 3: Adjacent residues reversed A/2 (GLU) then A/1 (MET) —
    // should be close with order_switch=true
    {
        std::pair<bool, bool> result =
            coot::peptide_C_and_N_are_close_p_gemmi(res2, res1);
        EXPECT_TRUE(result.first);   // are_close
        EXPECT_TRUE(result.second);  // order_switch
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}