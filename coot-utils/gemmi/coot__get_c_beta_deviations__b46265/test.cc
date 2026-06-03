#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find a residue by chain name and sequence number in a Structure
static const gemmi::Residue *find_residue(
    const gemmi::Structure &st, const std::string &chain_name, int seq_num)
{
    if (st.models.empty()) return nullptr;
    const gemmi::Model &model = st.models[0];
    for (const gemmi::Chain &chain : model.chains) {
        if (chain.name == chain_name) {
            for (const gemmi::Residue &res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, get_c_beta_deviations) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: full molecule — call get_c_beta_deviations on the Structure
    {
        std::map<const gemmi::Residue*, std::map<std::string, coot::c_beta_deviation_t>> results =
            coot::get_c_beta_deviations_gemmi(st);

        EXPECT_EQ(results.size(), 253u);

        // A/223 ASP
        {
            const gemmi::Residue *res = find_residue(st, "A", 223);
            ASSERT_NE(res, nullptr);
            EXPECT_STREQ(res->name.c_str(), "ASP");
            auto it = results.find(res);
            EXPECT_NE(it, results.end());
            if (it != results.end()) {
                auto const &alts = it->second;
                EXPECT_EQ(alts.size(), 1u);
                auto ait = alts.begin();
                EXPECT_STREQ(ait->first.c_str(), "");
                EXPECT_NEAR(ait->second.dist, 0.036081, 1e-4);
            }
        }

        // A/224 GLU
        {
            const gemmi::Residue *res = find_residue(st, "A", 224);
            ASSERT_NE(res, nullptr);
            EXPECT_STREQ(res->name.c_str(), "GLU");
            auto it = results.find(res);
            EXPECT_NE(it, results.end());
            if (it != results.end()) {
                auto const &alts = it->second;
                EXPECT_EQ(alts.size(), 1u);
                auto ait = alts.begin();
                EXPECT_STREQ(ait->first.c_str(), "");
                EXPECT_NEAR(ait->second.dist, 0.0573215, 1e-4);
            }
        }

        // A/225 VAL
        {
            const gemmi::Residue *res = find_residue(st, "A", 225);
            ASSERT_NE(res, nullptr);
            EXPECT_STREQ(res->name.c_str(), "VAL");
            auto it = results.find(res);
            EXPECT_NE(it, results.end());
            if (it != results.end()) {
                auto const &alts = it->second;
                EXPECT_EQ(alts.size(), 1u);
                auto ait = alts.begin();
                EXPECT_STREQ(ait->first.c_str(), "");
                EXPECT_NEAR(ait->second.dist, 0.0172193, 1e-4);
            }
        }

        // A/226 VAL
        {
            const gemmi::Residue *res = find_residue(st, "A", 226);
            ASSERT_NE(res, nullptr);
            EXPECT_STREQ(res->name.c_str(), "VAL");
            auto it = results.find(res);
            EXPECT_NE(it, results.end());
            if (it != results.end()) {
                auto const &alts = it->second;
                EXPECT_EQ(alts.size(), 1u);
                auto ait = alts.begin();
                EXPECT_STREQ(ait->first.c_str(), "");
                EXPECT_NEAR(ait->second.dist, 0.0756042, 1e-4);
            }
        }

        // A/227 TRP
        {
            const gemmi::Residue *res = find_residue(st, "A", 227);
            ASSERT_NE(res, nullptr);
            EXPECT_STREQ(res->name.c_str(), "TRP");
            auto it = results.find(res);
            EXPECT_NE(it, results.end());
            if (it != results.end()) {
                auto const &alts = it->second;
                EXPECT_EQ(alts.size(), 1u);
                auto ait = alts.begin();
                EXPECT_STREQ(ait->first.c_str(), "");
                EXPECT_NEAR(ait->second.dist, 0.00973383, 1e-4);
            }
        }
    }

    // Case 2: per-residue version for a specific residue with CB
    {
        gemmi::Residue *res = const_cast<gemmi::Residue*>(find_residue(st, "A", 10));
        ASSERT_NE(res, nullptr);
        std::map<std::string, coot::c_beta_deviation_t> res_results =
            coot::get_c_beta_deviations_residue_gemmi(*res);
        EXPECT_EQ(res_results.size(), 1u);
        for (auto const &alt : res_results) {
            EXPECT_STREQ(alt.first.c_str(), "");
            EXPECT_NEAR(alt.second.dist, 0.20076, 1e-4);
        }
    }

    // Case 3: edge case — residue //A/1 (MET)
    {
        gemmi::Residue *res = const_cast<gemmi::Residue*>(find_residue(st, "A", 1));
        ASSERT_NE(res, nullptr);
        EXPECT_STREQ(res->name.c_str(), "MET");
        std::map<std::string, coot::c_beta_deviation_t> res_results =
            coot::get_c_beta_deviations_residue_gemmi(*res);
        EXPECT_EQ(res_results.size(), 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}