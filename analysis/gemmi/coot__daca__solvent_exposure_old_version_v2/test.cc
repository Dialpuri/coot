#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, solvent_exposure_old_version_v2) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    {
        int n_models = static_cast<int>(st.models.size());
        EXPECT_EQ(n_models, 1);
    }

    // Case 1: side_chain_only = false
    {
        bool side_chain_only = false;
        auto result = coot::daca::solvent_exposure_old_version_v2_gemmi(&st, side_chain_only);
        EXPECT_EQ(result.size(), 267u);

        // First 5 residues
        EXPECT_NE(result[0].first.residue, nullptr);
        EXPECT_EQ(result[0].first.chain->name, "A");
        EXPECT_EQ(result[0].first.residue->seqid.num.value, 1);
        EXPECT_EQ(result[0].first.residue->name, "MET");
        EXPECT_NEAR(result[0].second, 43.0f, 1e-4);

        EXPECT_NE(result[1].first.residue, nullptr);
        EXPECT_EQ(result[1].first.chain->name, "A");
        EXPECT_EQ(result[1].first.residue->seqid.num.value, 2);
        EXPECT_EQ(result[1].first.residue->name, "GLU");
        EXPECT_NEAR(result[1].second, 28.0f, 1e-4);

        EXPECT_NE(result[2].first.residue, nullptr);
        EXPECT_EQ(result[2].first.chain->name, "A");
        EXPECT_EQ(result[2].first.residue->seqid.num.value, 3);
        EXPECT_EQ(result[2].first.residue->name, "ASN");
        EXPECT_NEAR(result[2].second, 52.0f, 1e-4);

        EXPECT_NE(result[3].first.residue, nullptr);
        EXPECT_EQ(result[3].first.chain->name, "A");
        EXPECT_EQ(result[3].first.residue->seqid.num.value, 4);
        EXPECT_EQ(result[3].first.residue->name, "PHE");
        EXPECT_NEAR(result[3].second, 67.0f, 1e-4);

        EXPECT_NE(result[4].first.residue, nullptr);
        EXPECT_EQ(result[4].first.chain->name, "A");
        EXPECT_EQ(result[4].first.residue->seqid.num.value, 5);
        EXPECT_EQ(result[4].first.residue->name, "GLN");
        EXPECT_NEAR(result[4].second, 49.0f, 1e-4);
    }

    // Case 2: side_chain_only = true
    {
        bool side_chain_only = true;
        auto result = coot::daca::solvent_exposure_old_version_v2_gemmi(&st, side_chain_only);
        EXPECT_EQ(result.size(), 253u);

        // First 5 residues
        EXPECT_NE(result[0].first.residue, nullptr);
        EXPECT_EQ(result[0].first.chain->name, "A");
        EXPECT_EQ(result[0].first.residue->seqid.num.value, 1);
        EXPECT_EQ(result[0].first.residue->name, "MET");
        EXPECT_NEAR(result[0].second, 28.0f, 1e-4);

        EXPECT_NE(result[1].first.residue, nullptr);
        EXPECT_EQ(result[1].first.chain->name, "A");
        EXPECT_EQ(result[1].first.residue->seqid.num.value, 2);
        EXPECT_EQ(result[1].first.residue->name, "GLU");
        EXPECT_NEAR(result[1].second, 17.0f, 1e-4);

        EXPECT_NE(result[2].first.residue, nullptr);
        EXPECT_EQ(result[2].first.chain->name, "A");
        EXPECT_EQ(result[2].first.residue->seqid.num.value, 3);
        EXPECT_EQ(result[2].first.residue->name, "ASN");
        EXPECT_NEAR(result[2].second, 29.0f, 1e-4);

        EXPECT_NE(result[3].first.residue, nullptr);
        EXPECT_EQ(result[3].first.chain->name, "A");
        EXPECT_EQ(result[3].first.residue->seqid.num.value, 4);
        EXPECT_EQ(result[3].first.residue->name, "PHE");
        EXPECT_NEAR(result[3].second, 52.0f, 1e-4);

        EXPECT_NE(result[4].first.residue, nullptr);
        EXPECT_EQ(result[4].first.chain->name, "A");
        EXPECT_EQ(result[4].first.residue->seqid.num.value, 5);
        EXPECT_EQ(result[4].first.residue->name, "GLN");
        EXPECT_NEAR(result[4].second, 29.0f, 1e-4);
    }

    // Case 3: null mol (edge case)
    {
        bool side_chain_only = false;
        auto result = coot::daca::solvent_exposure_old_version_v2_gemmi(nullptr, side_chain_only);
        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}