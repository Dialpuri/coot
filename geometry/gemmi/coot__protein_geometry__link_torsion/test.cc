#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"

TEST(GemmiLinkTorsionTest, two_rows) {
    // Build a gemmi::cif::Loop with two torsion rows
    gemmi::cif::Loop loop;
    loop.tags = {
        "link_id", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4",
        "atom_1_comp_id", "atom_2_comp_id", "atom_3_comp_id", "atom_4_comp_id",
        "value_angle", "value_angle_esd", "period"
    };

    // Row 0: phi torsion
    loop.add_row({"phi", "phi", "C", "N", "CA", "C", "1", "2", "2", "2", "-60", "5", "360"});

    // Row 1: psi torsion
    loop.add_row({"psi", "psi", "N", "CA", "C", "N", "3", "4", "4", "4", "-45", "3", "360"});

    EXPECT_EQ(loop.length(), 2u);

    EXPECT_NO_THROW(coot::protein_geometry::link_torsion_gemmi(loop));
}

TEST(GemmiLinkTorsionTest, empty_loop) {
    gemmi::cif::Loop empty_loop;
    empty_loop.tags = {"link_id"};
    // No rows added

    EXPECT_EQ(empty_loop.length(), 0u);

    EXPECT_NO_THROW(coot::protein_geometry::link_torsion_gemmi(empty_loop));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}