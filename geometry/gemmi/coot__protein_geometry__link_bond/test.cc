#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(OracleTest, link_bond_one_valid_row) {
    std::vector<coot::dictionary_residue_link_restraints_t> dict_link_res_restraints;

    gemmi::cif::Loop loop;
    loop.tags = {"link_id", "atom_id_1", "atom_id_2", "atom_1_comp_id", "atom_2_comp_id", "value_dist", "value_dist_esd"};
    loop.values = {"LINK1", "O", "C", "1", "2", "1.33", "0.01"};

    int n_bonds = coot::protein_geometry::link_bond_gemmi(loop, dict_link_res_restraints);

    EXPECT_EQ(loop.length(), 1);
    EXPECT_EQ(n_bonds, 1);
}

TEST(OracleTest, link_bond_empty_loop) {
    std::vector<coot::dictionary_residue_link_restraints_t> dict_link_res_restraints;

    gemmi::cif::Loop loop2;
    loop2.tags = {"link_id", "atom_id_1", "atom_id_2", "atom_1_comp_id", "atom_2_comp_id", "value_dist", "value_dist_esd"};
    // no values => length 0

    int n_bonds2 = coot::protein_geometry::link_bond_gemmi(loop2, dict_link_res_restraints);

    EXPECT_EQ(loop2.length(), 0);
    EXPECT_EQ(n_bonds2, 0);
}

TEST(OracleTest, link_bond_two_valid_rows) {
    std::vector<coot::dictionary_residue_link_restraints_t> dict_link_res_restraints;

    gemmi::cif::Loop loop3;
    loop3.tags = {"link_id", "atom_id_1", "atom_id_2", "atom_1_comp_id", "atom_2_comp_id", "value_dist", "value_dist_esd"};
    // row 0: LINK_A, O, N, 1, 3, 1.45, 0.02
    // row 1: LINK_B, C, CA, 2, 4, 1.52, 0.03
    loop3.values = {
        "LINK_A", "O", "N", "1", "3", "1.45", "0.02",
        "LINK_B", "C", "CA", "2", "4", "1.52", "0.03"
    };

    int n_bonds3 = coot::protein_geometry::link_bond_gemmi(loop3, dict_link_res_restraints);

    EXPECT_EQ(loop3.length(), 2);
    EXPECT_EQ(n_bonds3, 2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}