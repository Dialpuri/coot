#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "coot-utils/atom-tree.hh"
#include "coot-utils/atom-vertex.hh"
#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__construct_internal/gemmi/function.hh"

TEST(OracleTest, atom_tree_fill_torsions_gemmi) {
    // Read with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model &gemmi_model = st.models[0];
    gemmi::Chain *gemmi_chain = &gemmi_model.chains[0];
    gemmi::Residue *gemmi_res = &gemmi_chain->residues[0];

    ASSERT_NE(gemmi_res, nullptr);
    EXPECT_EQ(gemmi_chain->name, "A");
    EXPECT_EQ(gemmi_res->seqid.num.value, 1);
    EXPECT_EQ(gemmi_res->name, "MET");

    // Manually build restraints for MET
    coot::dictionary_residue_restraints_t rest("MET");

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

    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("N", "CA", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "C", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("C", "O", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CA", "CB", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CB", "CG", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("CG", "SD", ""));
    rest.bond_restraint.push_back(coot::dict_bond_restraint_t("SD", "CE", ""));

    rest.torsion_restraint.push_back(coot::dict_torsion_restraint_t("CHI1", "N", "CA", "CB", "CG", 60.0, 20.0, 3));
    rest.torsion_restraint.push_back(coot::dict_torsion_restraint_t("CHI2", "CA", "CB", "CG", "SD", 70.0, 20.0, 2));
    rest.torsion_restraint.push_back(coot::dict_torsion_restraint_t("CHI3", "CB", "CG", "SD", "CE", 90.0, 20.0, 2));

    EXPECT_EQ(rest.torsion_restraint.size(), 3u);
    EXPECT_EQ(rest.tree.size(), 8u);
    EXPECT_EQ(rest.bond_restraint.size(), 7u);

    // Build atom_vertex_vec using construct_internal_gemmi (no MMDB)
    std::vector<std::pair<int, int>> bonds_out;
    std::map<std::string, coot::map_index_t> name_to_index_out;
    std::vector<coot::atom_vertex> atom_vertex_vec_out;

    coot::construct_internal_gemmi(rest, gemmi_res, "", bonds_out, name_to_index_out, atom_vertex_vec_out);

    // Count torsions after construction (construct_internal calls fill_torsions internally
    // in the original MMDB code, but construct_internal_gemmi is gemmi-only and does NOT
    // call fill_torsions_gemmi — so we expect 0 torsions initially)
    int torsion_count = 0;
    for (size_t i = 0; i < atom_vertex_vec_out.size(); i++) {
        if (atom_vertex_vec_out[i].torsion_quad.first) {
            torsion_count++;
        }
    }

    EXPECT_EQ(atom_vertex_vec_out.size(), 8u);
    EXPECT_EQ(torsion_count, 0);

    // Call fill_torsions_gemmi directly with gemmi residue
    bool result = coot::fill_torsions_gemmi(atom_vertex_vec_out, rest, gemmi_res, "");
    EXPECT_TRUE(result);

    // Count again after direct call
    int torsion_count_after = 0;
    for (size_t i = 0; i < atom_vertex_vec_out.size(); i++) {
        if (atom_vertex_vec_out[i].torsion_quad.first) {
            torsion_count_after++;
        }
    }
    EXPECT_EQ(torsion_count_after, 3);

    // Case 2: empty torsion restraints - should return false
    coot::dictionary_residue_restraints_t empty_rest;
    bool result2 = coot::fill_torsions_gemmi(atom_vertex_vec_out, empty_rest, gemmi_res, "");
    EXPECT_FALSE(result2);

    // Case 3: torsions referencing non-existent atoms - should return false
    {
        coot::dictionary_residue_restraints_t bad_rest("MET");
        bad_rest.tree.push_back(t0); // copy the START node
        bad_rest.torsion_restraint.push_back(coot::dict_torsion_restraint_t("BAD", "NON1", "NON2", "NON3", "NON4", 0.0, 0.0, 1));

        std::vector<std::pair<int, int>> bad_bonds_out;
        std::map<std::string, coot::map_index_t> bad_name_to_index_out;
        std::vector<coot::atom_vertex> bad_atom_vertex_vec_out;

        coot::construct_internal_gemmi(bad_rest, gemmi_res, "", bad_bonds_out, bad_name_to_index_out, bad_atom_vertex_vec_out);

        bool result3 = coot::fill_torsions_gemmi(bad_atom_vertex_vec_out, bad_rest, gemmi_res, "");
        EXPECT_FALSE(result3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
