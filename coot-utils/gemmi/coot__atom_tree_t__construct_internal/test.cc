#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "coot-utils/atom-tree.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (auto& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (auto& res : chain.residues) {
            if (res.seqid.num.value == seq_num) {
                return &res;
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, atom_tree_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Use the first model
    gemmi::Model& model = st.models[0];

    // --- Case 1: Valid residue (MET 1, chain A) ---
    {
        gemmi::Residue* res = find_residue(model, "A", 1);
        EXPECT_NE(res, nullptr);
        EXPECT_EQ(res->seqid.num.value, 1);
        EXPECT_EQ(res->name, "MET");

        coot::dictionary_residue_restraints_t rest("GLY");

        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("N", "CA", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "C", ""));
        rest.bond_restraint.push_back(coot::dict_bond_restraint_t("C", "O", ""));

        coot::dict_chem_comp_tree_t t0("N", "", "", "");
        t0.connect_type = "START";
        rest.tree.push_back(t0);
        rest.tree.push_back(coot::dict_chem_comp_tree_t("CA", "N", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("C", "CA", "", ""));
        rest.tree.push_back(coot::dict_chem_comp_tree_t("O", "C", "", ""));

        EXPECT_EQ(rest.bond_restraint.size(), 3u);
        EXPECT_EQ(rest.tree.size(), 4u);

        std::vector<std::pair<int,int>> bonds;
        std::map<std::string, coot::map_index_t> name_to_index;
        std::vector<coot::atom_vertex> atom_vertex_vec;

        coot::construct_internal_gemmi(rest, res, "", bonds, name_to_index, atom_vertex_vec);

        EXPECT_EQ(bonds.size(), 3u);
        EXPECT_EQ(atom_vertex_vec.size(), 19u);
        EXPECT_EQ(name_to_index.size(), 20u);
    }

    // --- Case 2: Non-existent residue (should fail to find residue) ---
    {
        gemmi::Residue* res = find_residue(model, "Z", 999);
        EXPECT_EQ(res, nullptr);
    }

    // --- Case 3: Another residue with valid tree (GLU 2 in chain A) ---
    {
        gemmi::Residue* res = find_residue(model, "A", 2);
        if (!res) {
            EXPECT_FALSE(true); // should not happen based on oracle
        } else {
            EXPECT_EQ(res->seqid.num.value, 2);
            EXPECT_EQ(res->name, "GLU");

            coot::dictionary_residue_restraints_t rest("ALA");

            rest.bond_restraint.push_back(coot::dict_bond_restraint_t("N", "CA", ""));
            rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "C", ""));
            rest.bond_restraint.push_back(coot::dict_bond_restraint_t("C", "O", ""));
            rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "CB", ""));

            rest.tree.push_back(coot::dict_chem_comp_tree_t("N", "", "", ""));
            rest.tree.push_back(coot::dict_chem_comp_tree_t("CA", "N", "", ""));
            rest.tree.push_back(coot::dict_chem_comp_tree_t("C", "CA", "", ""));
            rest.tree.push_back(coot::dict_chem_comp_tree_t("O", "C", "", ""));
            rest.tree.push_back(coot::dict_chem_comp_tree_t("CB", "CA", "", ""));
            // Must set START on root atom for fill_atom_vertex_vec to work
            rest.tree[0].connect_type = "START";

            EXPECT_EQ(rest.bond_restraint.size(), 4u);

            std::vector<std::pair<int,int>> bonds;
            std::map<std::string, coot::map_index_t> name_to_index;
            std::vector<coot::atom_vertex> atom_vertex_vec;

            coot::construct_internal_gemmi(rest, res, "", bonds, name_to_index, atom_vertex_vec);

            EXPECT_EQ(bonds.size(), 4u);
            EXPECT_EQ(atom_vertex_vec.size(), 15u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}