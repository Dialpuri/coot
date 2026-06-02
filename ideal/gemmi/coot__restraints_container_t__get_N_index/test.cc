#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_N_index) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    auto& model = st.models[0];

    // Find chain A
    gemmi::Chain* chain_a = nullptr;
    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            chain_a = &chain;
            break;
        }
    }
    ASSERT_NE(chain_a, nullptr);

    // Case 1: valid residue with N atom — residue 1 of chain A (MET)
    {
        gemmi::Residue* res = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 1) {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        EXPECT_EQ(res->name, "MET");
        EXPECT_EQ(chain_a->name, "A");
        EXPECT_EQ(res->seqid.num.value, 1);

        int n_index = coot::restraints_container_t::get_N_index_gemmi(*res);
        EXPECT_EQ(n_index, 0);
    }

    // Case 2: another residue — residue 10 of chain A (ILE)
    {
        gemmi::Residue* res2 = nullptr;
        for (auto& r : chain_a->residues) {
            if (r.seqid.num.value == 10) {
                res2 = &r;
                break;
            }
        }
        ASSERT_NE(res2, nullptr);

        EXPECT_EQ(res2->name, "ILE");
        EXPECT_EQ(chain_a->name, "A");
        EXPECT_EQ(res2->seqid.num.value, 10);

        int n_index = coot::restraints_container_t::get_N_index_gemmi(*res2);
        EXPECT_EQ(n_index, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}