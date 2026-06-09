#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t) {
    // Case 1: Typical amino acid residue
    {
        std::string comp_id = "ALA";
        int read_number = 1;
        coot::dictionary_residue_restraints_t rest =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(comp_id, read_number);

        EXPECT_EQ(rest.residue_info.comp_id, "ALA");
        EXPECT_EQ(rest.read_number, 1);
        EXPECT_EQ(rest.has_partial_charges_flag, 0);
        EXPECT_EQ(rest.filled_with_bond_order_data_only_flag, 0);
        EXPECT_EQ(rest.nuclear_distances_flag, 0);
        EXPECT_EQ(rest.comp_id(), "ALA");
        EXPECT_EQ(rest.is_filled(), 0);
        EXPECT_EQ(rest.number_of_atoms(), 0u);
    }

    // Case 2: Ligand with different comp_id and read_number
    {
        std::string comp_id = "ATP";
        int read_number = 42;
        coot::dictionary_residue_restraints_t rest =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(comp_id, read_number);

        EXPECT_EQ(rest.residue_info.comp_id, "ATP");
        EXPECT_EQ(rest.read_number, 42);
        EXPECT_EQ(rest.has_partial_charges_flag, 0);
        EXPECT_EQ(rest.filled_with_bond_order_data_only_flag, 0);
        EXPECT_EQ(rest.nuclear_distances_flag, 0);
        EXPECT_EQ(rest.comp_id(), "ATP");
        EXPECT_EQ(rest.is_filled(), 0);
        EXPECT_EQ(rest.number_of_atoms(), 0u);
    }

    // Case 3: Empty string comp_id and negative read_number (edge case)
    {
        std::string comp_id = "";
        int read_number = -1;
        coot::dictionary_residue_restraints_t rest =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(comp_id, read_number);

        EXPECT_EQ(rest.residue_info.comp_id, "");
        EXPECT_EQ(rest.read_number, -1);
        EXPECT_EQ(rest.has_partial_charges_flag, 0);
        EXPECT_EQ(rest.filled_with_bond_order_data_only_flag, 0);
        EXPECT_EQ(rest.nuclear_distances_flag, 0);
        EXPECT_EQ(rest.comp_id(), "");
        EXPECT_EQ(rest.is_filled(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
