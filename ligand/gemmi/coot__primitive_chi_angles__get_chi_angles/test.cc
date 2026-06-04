#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number in model 0
static gemmi::Residue* find_residue(gemmi::Structure& st,
                                    const std::string& chain_name,
                                    int seq_num) {
    if (st.models.empty()) return nullptr;
    gemmi::Model& model = st.models[0];
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_name) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTestGemmi, primitive_chi_angles) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: PHE 4 — has chi angles (chi1, chi2)
    {
        gemmi::Residue* res = find_residue(st, "A", 4);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "PHE");
        EXPECT_EQ(res->seqid.num.value, 4);

        std::vector<coot::alt_confed_chi_angles> chis =
            coot::primitive_chi_angles::get_chi_angles_gemmi(*res);

        EXPECT_EQ(chis.size(), 1u);
        EXPECT_EQ(chis[0].alt_conf, "");
        EXPECT_EQ(chis[0].chi_angles.size(), 2u);
        EXPECT_EQ(chis[0].chi_angles[0].first, 1);
        EXPECT_NEAR(chis[0].chi_angles[0].second, -56.5473, 1e-4);
        EXPECT_EQ(chis[0].chi_angles[1].first, 2);
        EXPECT_NEAR(chis[0].chi_angles[1].second, -91.3803, 1e-4);
    }

    // Case 2: GLY 11 — no chi angles, should throw
    {
        gemmi::Residue* res = find_residue(st, "A", 11);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "GLY");
        EXPECT_EQ(res->seqid.num.value, 11);

        EXPECT_THROW({
            coot::primitive_chi_angles::get_chi_angles_gemmi(*res);
        }, std::runtime_error);
    }

    // Case 3: VAL 7 — has one chi angle (chi1)
    {
        gemmi::Residue* res = find_residue(st, "A", 7);
        ASSERT_NE(res, nullptr);
        EXPECT_EQ(res->name, "VAL");
        EXPECT_EQ(res->seqid.num.value, 7);

        std::vector<coot::alt_confed_chi_angles> chis =
            coot::primitive_chi_angles::get_chi_angles_gemmi(*res);

        EXPECT_EQ(chis.size(), 1u);
        EXPECT_EQ(chis[0].alt_conf, "");
        EXPECT_EQ(chis[0].chi_angles.size(), 1u);
        EXPECT_EQ(chis[0].chi_angles[0].first, 1);
        EXPECT_NEAR(chis[0].chi_angles[0].second, 177.07, 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}