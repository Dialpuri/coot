#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_daca_fragments) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_LT(0u, st.models.size());
    gemmi::Model& model = st.models[0];

    ASSERT_LT(0u, model.chains.size());
    gemmi::Chain& chain = model.chains[0];

    // Case 1: first residue (MET)
    {
        ASSERT_LT(0u, chain.residues.size());
        gemmi::Residue& res = chain.residues[0];
        EXPECT_EQ(res.name, "MET");
        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(res.seqid.num.value, 1);

        std::vector<std::vector<const gemmi::Atom*>> fragments =
            coot::get_daca_fragments_gemmi(&res);

        EXPECT_EQ(fragments.size(), 5u);

        // fragment 0: 3 atoms — CA, C, O
        EXPECT_EQ(fragments[0].size(), 3u);
        EXPECT_EQ(fragments[0][0]->name, "CA");
        EXPECT_EQ(fragments[0][1]->name, "C");
        EXPECT_EQ(fragments[0][2]->name, "O");

        // fragment 1: 4 atoms — N, C, CA, CB
        EXPECT_EQ(fragments[1].size(), 4u);
        EXPECT_EQ(fragments[1][0]->name, "N");
        EXPECT_EQ(fragments[1][1]->name, "C");
        EXPECT_EQ(fragments[1][2]->name, "CA");
        EXPECT_EQ(fragments[1][3]->name, "CB");

        // fragment 2: 3 atoms — CA, CB, CG
        EXPECT_EQ(fragments[2].size(), 3u);
        EXPECT_EQ(fragments[2][0]->name, "CA");
        EXPECT_EQ(fragments[2][1]->name, "CB");
        EXPECT_EQ(fragments[2][2]->name, "CG");

        // fragment 3: 3 atoms — CB, CG, SD
        EXPECT_EQ(fragments[3].size(), 3u);
        EXPECT_EQ(fragments[3][0]->name, "CB");
        EXPECT_EQ(fragments[3][1]->name, "CG");
        EXPECT_EQ(fragments[3][2]->name, "SD");

        // fragment 4: 3 atoms — CG, SD, CE
        EXPECT_EQ(fragments[4].size(), 3u);
        EXPECT_EQ(fragments[4][0]->name, "CG");
        EXPECT_EQ(fragments[4][1]->name, "SD");
        EXPECT_EQ(fragments[4][2]->name, "CE");
    }

    // Case 2: another residue (GLU)
    {
        if (chain.residues.size() > 1) {
            gemmi::Residue& res2 = chain.residues[1];
            EXPECT_EQ(res2.name, "GLU");
            EXPECT_EQ(chain.name, "A");
            EXPECT_EQ(res2.seqid.num.value, 2);

            std::vector<std::vector<const gemmi::Atom*>> fragments2 =
                coot::get_daca_fragments_gemmi(&res2);
            EXPECT_EQ(fragments2.size(), 5u);
        }
    }

    // Case 3: residue from chain B (HOH)
    {
        if (model.chains.size() > 1) {
            gemmi::Chain& chain2 = model.chains[1];
            if (!chain2.residues.empty()) {
                gemmi::Residue& res3 = chain2.residues[0];
                EXPECT_EQ(res3.name, "HOH");
                EXPECT_EQ(chain2.name, "B");
                EXPECT_EQ(res3.seqid.num.value, 1);

                std::vector<std::vector<const gemmi::Atom*>> fragments3 =
                    coot::get_daca_fragments_gemmi(&res3);
                EXPECT_EQ(fragments3.size(), 0u);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}