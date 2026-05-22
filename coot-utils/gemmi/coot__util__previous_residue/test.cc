#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, previous_residue_function) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // case 1: middle residue with a valid previous residue
    {
        gemmi::Chain* chain_p = nullptr;
        gemmi::Residue* target = nullptr;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 10) {
                    chain_p = &chain;
                    target = &res;
                    break;
                }
            }
        }
        ASSERT_NE(chain_p, nullptr);
        ASSERT_NE(target, nullptr);

        gemmi::CRA cra{chain_p, target, nullptr};
        gemmi::Residue* prev = coot::util::previous_residue_gemmi(cra);
        EXPECT_NE(prev, nullptr);
        if (prev) {
            EXPECT_EQ(prev->seqid.num.value, 9);
        }
    }

    // case 2: first residue — guard clause triggers (no previous)
    {
        gemmi::Chain* chain_p = nullptr;
        gemmi::Residue* target = nullptr;
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == 1) {
                    chain_p = &chain;
                    target = &res;
                    break;
                }
            }
        }
        ASSERT_NE(chain_p, nullptr);
        ASSERT_NE(target, nullptr);

        gemmi::CRA cra{chain_p, target, nullptr};
        gemmi::Residue* prev = coot::util::previous_residue_gemmi(cra);
        EXPECT_EQ(prev, nullptr);
    }

    // case 3: null input
    {
        gemmi::CRA cra{nullptr, nullptr, nullptr};
        gemmi::Residue* prev = coot::util::previous_residue_gemmi(cra);
        EXPECT_EQ(prev, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}