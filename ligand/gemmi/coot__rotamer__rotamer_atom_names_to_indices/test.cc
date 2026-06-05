#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, rotamer_atom_names_to_indices) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    // Case 1: GLY at //A/16 — no sidechain chi angles (should find 0 chi atom indices)
    {
        gemmi::Residue* gly_res = nullptr;
        gemmi::Chain* gly_chain = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                gly_chain = &chain;
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 16 && res.name == "GLY") {
                        gly_res = &res;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(gly_res, nullptr);

        EXPECT_EQ(gly_res->name, "GLY");
        EXPECT_EQ(gly_res->atoms.size(), 4u);

        std::vector<std::vector<std::string>> rotamer_ats;
        EXPECT_EQ(rotamer_ats.size(), 0u);

        std::vector<std::vector<int>> atom_indices =
            coot::rotamer::rotamer_atom_names_to_indices_gemmi(rotamer_ats, {gly_chain, gly_res, nullptr});

        EXPECT_EQ(atom_indices.size(), 0u);
    }

    // Case 2: VAL at //A/7 — has sidechain (should return 1 matched chi set)
    {
        gemmi::Residue* val_res = nullptr;
        gemmi::Chain* val_chain = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                val_chain = &chain;
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == 7 && res.name == "VAL") {
                        val_res = &res;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(val_res, nullptr);

        EXPECT_EQ(val_res->name, "VAL");
        EXPECT_EQ(val_res->atoms.size(), 7u);

        // VAL chi1: N, CA, CB, CG1
        std::vector<std::vector<std::string>> rotamer_ats;
        rotamer_ats.push_back({"N", "CA", "CB", "CG1"});
        EXPECT_EQ(rotamer_ats.size(), 1u);

        std::vector<std::vector<int>> atom_indices =
            coot::rotamer::rotamer_atom_names_to_indices_gemmi(rotamer_ats, {val_chain, val_res, nullptr});

        EXPECT_EQ(atom_indices.size(), 1u);
    }

    // Case 3: n_residue_atoms == 0 (edge case — should return empty)
    {
        std::vector<std::vector<std::string>> rotamer_ats;
        rotamer_ats.push_back({"N", "CA", "C", "O"});

        EXPECT_EQ(rotamer_ats.size(), 1u);

        std::vector<std::vector<int>> atom_indices =
            coot::rotamer::rotamer_atom_names_to_indices_gemmi(rotamer_ats, {nullptr, nullptr, nullptr});

        EXPECT_EQ(atom_indices.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}