#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTestGemmi, atom_tree_t) {
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // Find chain A, residue 1 (MET)
    gemmi::Residue *res = nullptr;
    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != "A") continue;
        for (gemmi::Residue &r : chain.residues) {
            if (r.seqid.num.value == 1 && r.name == "MET") {
                res = &r;
                break;
            }
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr);

    // Verify the residue properties
    ASSERT_EQ(res->seqid.num.value, 1);
    ASSERT_EQ(res->name, "MET");

    // Case 1: Normal call
    {
        coot::dictionary_residue_restraints_t rest("MET");

        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("N", "CA", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "C", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("C", "O", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "CB", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CB", "CG", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CG", "SD", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("SD", "CE", ""));

        coot::dict_chem_comp_tree_t t0("N", "", "", "");
        t0.connect_type = "START";
        rest.tree.push_back(t0);
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CA", "N", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("C", "CA", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("O", "C", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CB", "CA", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CG", "CB", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("SD", "CG", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CE", "SD", "", ""));

        EXPECT_EQ(rest.comp_id(), "");
        EXPECT_EQ(rest.tree.size(), 8u);
        EXPECT_EQ(rest.bond_restraint.size(), 7u);

        // Call construct_internal_gemmi directly
        std::vector<std::pair<int, int>> bonds;
        std::map<std::string, coot::map_index_t> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        coot::construct_internal_gemmi(rest, res, "", bonds, name_to_index, atom_vertex_vec);

        EXPECT_EQ(bonds.size(), 7u);
        EXPECT_EQ(atom_vertex_vec.size(), 19u);
        EXPECT_EQ(name_to_index.size(), 20u);
    }

    // Case 2: Null residue — guard clause should throw
    {
        coot::dictionary_residue_restraints_t rest("MET");
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("N", "CA", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("N", "", "", ""));
        rest.tree[0].connect_type = "START";
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CA", "N", "", ""));

        std::vector<std::pair<int, int>> bonds;
        std::map<std::string, coot::map_index_t> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        EXPECT_THROW(
            coot::construct_internal_gemmi(rest, nullptr, "", bonds, name_to_index, atom_vertex_vec),
            std::runtime_error
        );
    }

    // Case 3: Empty restraints (no tree) — should throw
    {
        coot::dictionary_residue_restraints_t empty_rest;
        EXPECT_EQ(empty_rest.comp_id(), "");
        EXPECT_EQ(empty_rest.tree.size(), 0u);

        std::vector<std::pair<int, int>> bonds;
        std::map<std::string, coot::map_index_t> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        EXPECT_THROW(
            coot::construct_internal_gemmi(empty_rest, res, "", bonds, name_to_index, atom_vertex_vec),
            std::runtime_error
        );
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
