#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t) {
    // Case 1: default constructor
    {
        auto rest = coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi();
        EXPECT_EQ(rest.filled_with_bond_order_data_only_flag, 0);
        EXPECT_EQ(rest.has_partial_charges_flag, 0);
        EXPECT_EQ(rest.read_number, -1);
        EXPECT_EQ(rest.nuclear_distances_flag, 0);
    }

    // Case 2: (const std::string &comp_id_in, int read_number_in)
    {
        auto rest = coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi("ALA", 42);
        EXPECT_EQ(rest.residue_info.comp_id, "ALA");
        EXPECT_EQ(rest.read_number, 42);
        EXPECT_EQ(rest.has_partial_charges_flag, 0);
        EXPECT_EQ(rest.filled_with_bond_order_data_only_flag, 0);
        EXPECT_EQ(rest.nuclear_distances_flag, 0);
    }

    // Case 3: (bool constructor_for_srs_restraints) — 7 converts to true
    {
        auto rest = coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(7);
        EXPECT_EQ(rest.read_number, -1);
        EXPECT_EQ(rest.residue_info.comp_id, "");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
