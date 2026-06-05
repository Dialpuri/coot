#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, O3prime_and_P_are_close_ng) {
    // Read the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model &model = st.models[0];

    // Find chain A
    const gemmi::Chain *chain_a = nullptr;
    for (auto &ch : model.chains) {
        if (ch.name == "A") {
            chain_a = &ch;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Find residue 10 and 11 from chain A
    const gemmi::Residue *res1 = nullptr;
    const gemmi::Residue *res2 = nullptr;
    for (auto &res : chain_a->residues) {
        if (res.seqid.num.value == 10) {
            res1 = &res;
        }
        if (res.seqid.num.value == 11) {
            res2 = &res;
        }
    }

    ASSERT_NE(res1, nullptr);
    ASSERT_NE(res2, nullptr);

    // Case 1: Two protein residues — O3' and P atoms don't exist → false
    {
        float d_crit = 2.0f;
        bool result = coot::O3prime_and_P_are_close_ng_gemmi(
            const_cast<gemmi::Residue*>(res1),
            const_cast<gemmi::Residue*>(res2),
            d_crit);
        EXPECT_EQ(result, false);
    }

    // Case 2: Same residues, very large d_crit — still no O3'/P atoms → false
    {
        float d_crit = 100.0f;
        bool result = coot::O3prime_and_P_are_close_ng_gemmi(
            const_cast<gemmi::Residue*>(res1),
            const_cast<gemmi::Residue*>(res2),
            d_crit);
        EXPECT_EQ(result, false);
    }

    // Case 3: Check that O3' and P atoms are actually absent
    {
        bool o3_found = false;
        for (const auto &atom : res1->atoms) {
            if (atom.name == " O3'") {
                o3_found = true;
                break;
            }
        }
        EXPECT_EQ(o3_found, false);

        bool p_found = false;
        for (const auto &atom : res2->atoms) {
            if (atom.name == " P  ") {
                p_found = true;
                break;
            }
        }
        EXPECT_EQ(p_found, false);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}