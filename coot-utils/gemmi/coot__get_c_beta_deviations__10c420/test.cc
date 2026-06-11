#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Trim whitespace helper
static std::string trim_str(const std::string &s) {
   size_t start = s.find_first_not_of(" \t\r\n");
   if (start == std::string::npos) return "";
   size_t end = s.find_last_not_of(" \t\r\n");
   return s.substr(start, end - start + 1);
}

TEST(OracleTest, get_c_beta_deviations) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model *model_p = &(st.models[0]);
    ASSERT_NE(model_p, nullptr);

    gemmi::Chain *chain_p = &(model_p->chains[0]);
    ASSERT_NE(chain_p, nullptr);

    // Case from oracle: Glycine (no CB) — should return empty map
    {
        gemmi::Residue *gly_res = nullptr;
        for (gemmi::Residue &res : chain_p->residues) {
            if (trim_str(res.name) == "GLY") {
                gly_res = &res;
                break;
            }
        }
        ASSERT_NE(gly_res, nullptr);
        EXPECT_EQ(std::string(chain_p->name), std::string("A"));
        EXPECT_EQ(gly_res->seqid.num.value, 11);

        std::map<std::string, coot::c_beta_deviation_gemmi> result = coot::get_c_beta_deviations_gemmi(gly_res);
        EXPECT_EQ(result.size(), 0u);
    }

    // Complementary case: MET at seqnum 1 has CB atom — should produce non-empty result
    {
        gemmi::Residue *met_res = nullptr;
        for (gemmi::Residue &res : chain_p->residues) {
            if (trim_str(res.name) == "MET") {
                met_res = &res;
                break;
            }
        }
        ASSERT_NE(met_res, nullptr);
        EXPECT_EQ(met_res->seqid.num.value, 1);

        std::map<std::string, coot::c_beta_deviation_gemmi> result = coot::get_c_beta_deviations_gemmi(met_res);
        EXPECT_GT(result.size(), 0u);

        // Check that the main (empty string) alt conf has meaningful data
        auto it = result.find("");
        EXPECT_NE(it, result.end());
        EXPECT_NE(it->second.at, nullptr);
        EXPECT_FALSE(it->second.pos_ideal.is_null());
        EXPECT_GE(it->second.dist, 0.0);
    }

    // Case 3 from oracle: no partial residue found (all residues have >= 3 atoms)
    {
        gemmi::Residue *partial_res = nullptr;
        for (gemmi::Residue &res : chain_p->residues) {
            if (res.atoms.size() < 3) {
                partial_res = &res;
                break;
            }
        }
        EXPECT_EQ(partial_res, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}