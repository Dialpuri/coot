#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, calc_cablam) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model &model = st.models[0];
    ASSERT_FALSE(model.chains.empty());

    gemmi::Chain &chain = model.chains[0];

    int n_res = static_cast<int>(chain.residues.size());
    EXPECT_EQ(n_res, 267);

    // Case 1: typical interior residue (ires=10)
    {
        int ires = 10;
        gemmi::Residue &residue_this = chain.residues[ires];
        EXPECT_EQ(chain.name, "A");
        EXPECT_EQ(ires, 10);
        EXPECT_EQ(residue_this.name, "GLY");
        EXPECT_EQ(residue_this.seqid.num.value, 11);

        coot::cablam_markup_t_gemmi cm = coot::calc_cablam_gemmi(&chain, &residue_this, ires, 0.5);
        ASSERT_NE(cm.residue, nullptr);
        EXPECT_EQ(cm.residue->name, "GLY");
        EXPECT_NEAR(cm.score, 0.5, 1e-4);
    }

    // Case 2: near start (ires=1)
    {
        int ires = 1;
        gemmi::Residue &residue_this = chain.residues[ires];
        EXPECT_EQ(ires, 1);
        EXPECT_EQ(residue_this.name, "GLU");
        EXPECT_EQ(residue_this.seqid.num.value, 2);

        coot::cablam_markup_t_gemmi cm = coot::calc_cablam_gemmi(&chain, &residue_this, ires, 0.3);
        ASSERT_NE(cm.residue, nullptr);
        EXPECT_EQ(cm.residue->name, "GLU");
        EXPECT_NEAR(cm.score, 0.3, 1e-4);
    }

    // Case 3: near end (ires=n_res-3=264)
    {
        int ires = n_res - 3;
        gemmi::Residue &residue_this = chain.residues[ires];
        EXPECT_EQ(ires, 264);
        EXPECT_EQ(residue_this.name, "LEU");
        EXPECT_EQ(residue_this.seqid.num.value, 296);

        coot::cablam_markup_t_gemmi cm = coot::calc_cablam_gemmi(&chain, &residue_this, ires, 0.7);
        ASSERT_NE(cm.residue, nullptr);
        EXPECT_EQ(cm.residue->name, "LEU");
        EXPECT_NEAR(cm.score, 0.7, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}