#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t) {
    // Case 1: constructor with false
    {
        bool constructor_for_srs = false;
        coot::dictionary_residue_restraints_t rest =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(constructor_for_srs);

        EXPECT_TRUE(rest.filled_with_bond_order_data_only_flag);
        EXPECT_FALSE(rest.has_partial_charges_flag);
        EXPECT_EQ(rest.read_number, -1);
        EXPECT_FALSE(rest.nuclear_distances_flag);
    }

    // Case 2: constructor with true
    {
        bool constructor_for_srs = true;
        coot::dictionary_residue_restraints_t rest =
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(constructor_for_srs);

        EXPECT_TRUE(rest.filled_with_bond_order_data_only_flag);
        EXPECT_FALSE(rest.has_partial_charges_flag);
        EXPECT_EQ(rest.read_number, -1);
        EXPECT_FALSE(rest.nuclear_distances_flag);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
