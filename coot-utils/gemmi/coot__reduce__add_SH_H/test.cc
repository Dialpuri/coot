#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_SH_H) {
    {
        // Case 1: CYS at //A/7 — has SG atom, should add 1 H
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/canyon/test/canyon-test.pdb");
        gemmi::Model &model = st.models[0];

        // Find chain A, residue 7
        gemmi::Chain *chain_p = nullptr;
        gemmi::Residue *res_cys = nullptr;
        for (auto &chain : model.chains) {
            if (chain.name == "A") {
                chain_p = &chain;
                for (auto &res : chain.residues) {
                    if (res.seqid.num.value == 7) {
                        res_cys = &res;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res_cys, nullptr);

        int atoms_before = static_cast<int>(res_cys->atoms.size());
        EXPECT_EQ(atoms_before, 6);

        gemmi::CRA cra{chain_p, res_cys, nullptr};

        std::vector<gemmi::Atom *> result = coot::reduce::add_SH_H_gemmi(
            " HG ", " SG ", " CB ", " CA ", 0.2, 109.5, 180, *res_cys, cra, model, st);

        EXPECT_EQ(result.size(), 1);
        EXPECT_EQ(res_cys->atoms.size(), 7);
        ASSERT_FALSE(result.empty());
        ASSERT_NE(result[0], nullptr);
        EXPECT_EQ(std::string(result[0]->name), " HG ");
    }

    {
        // Case 2: ASP at //A/1 — no SG atom, should return empty
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/canyon/test/canyon-test.pdb");
        gemmi::Model &model = st.models[0];

        // Find chain A, residue 1
        gemmi::Chain *chain_p = nullptr;
        gemmi::Residue *res_ala = nullptr;
        for (auto &chain : model.chains) {
            if (chain.name == "A") {
                chain_p = &chain;
                for (auto &res : chain.residues) {
                    if (res.seqid.num.value == 1) {
                        res_ala = &res;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res_ala, nullptr);

        int atoms_before = static_cast<int>(res_ala->atoms.size());
        EXPECT_EQ(atoms_before, 8);

        gemmi::CRA cra{chain_p, res_ala, nullptr};

        std::vector<gemmi::Atom *> result = coot::reduce::add_SH_H_gemmi(
            " HG ", " SG ", " CB ", " CA ", 0.2, 109.5, 180, *res_ala, cra, model, st);

        EXPECT_EQ(result.size(), 0);
        EXPECT_EQ(res_ala->atoms.size(), 8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}