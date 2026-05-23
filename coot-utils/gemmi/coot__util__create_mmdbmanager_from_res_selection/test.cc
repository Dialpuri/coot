#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, create_mmdbmanager_from_res_selection) {
    // Load the CIF file using gemmi (auto-detects format)
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty()) << "Failed to load structure";

    // Case 1: Single residue (ligand LZA 1299)
    {
        std::vector<gemmi::Residue*> sel_res;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 1299) {
                        sel_res.push_back(&res);
                        break;
                    }
                }
            }
        }
        EXPECT_EQ(sel_res.size(), 1u) << "Could not find residue 1299";

        auto result = coot::util::create_mmdbmanager_from_res_selection_gemmi(
            st, sel_res, 0, 0, "", "A", 0
        );

        EXPECT_FALSE(result.first.models.empty()) << "OUTPUT new_mol: created";
        EXPECT_EQ(result.second, 16777217) << "OUTPUT atom_index_handle";
    }

    // Case 2: Multiple consecutive residues (10, 11, 12)
    {
        std::vector<gemmi::Residue*> sel_res;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10 || res.seqid.num.value == 11 || res.seqid.num.value == 12) {
                        sel_res.push_back(&res);
                    }
                }
            }
        }
        EXPECT_GE(sel_res.size(), 3u) << "Could not find residues 10, 11, 12";

        auto result = coot::util::create_mmdbmanager_from_res_selection_gemmi(
            st, sel_res, 0, 0, "", "A", 0
        );

        EXPECT_EQ(result.first.models[0].chains[0].residues.size(), 3u) << "OUTPUT new_mol.n_residues";
    }

    // Case 3: With altconf specification
    {
        std::vector<gemmi::Residue*> sel_res;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        sel_res.push_back(&res);
                        break;
                    }
                }
            }
        }
        EXPECT_EQ(sel_res.size(), 1u) << "Could not find residue 10";

        auto result = coot::util::create_mmdbmanager_from_res_selection_gemmi(
            st, sel_res, 0, 0, "A", "A", 0
        );

        EXPECT_FALSE(result.first.models.empty()) << "OUTPUT new_mol: created";
    }

    // Case 4: With whole_res_flag (residue_from_alt_conf_split_flag = 1)
    {
        std::vector<gemmi::Residue*> sel_res;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        sel_res.push_back(&res);
                        break;
                    }
                }
            }
        }
        EXPECT_EQ(sel_res.size(), 1u) << "Could not find residue 10";

        auto result = coot::util::create_mmdbmanager_from_res_selection_gemmi(
            st, sel_res, 0, 0, "", "A", 1
        );

        EXPECT_FALSE(result.first.models.empty()) << "OUTPUT new_mol: created";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}