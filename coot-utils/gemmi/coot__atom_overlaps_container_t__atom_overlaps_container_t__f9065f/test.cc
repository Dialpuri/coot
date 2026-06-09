#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps) {
    // Case 1: default constructor — the target function
    {
        coot::atom_overlaps_container_t ao = coot::atom_overlaps_container_t::atom_overlaps_container_t_gemmi();
        EXPECT_EQ(ao.overlaps.size(), 0u);
        EXPECT_FALSE(ao.get_have_dictionary());
    }

    // Case 2: atom_overlap_t default constructor (nested type, same header)
    {
        coot::atom_overlap_t aot;
        EXPECT_NEAR(aot.overlap_volume, -1.0, 1e-4);
        EXPECT_FALSE(aot.is_h_bond);
        EXPECT_EQ(aot.ligand_atom_index, -1);
    }

    // Case 3: atom_overlaps_dots_container_t default constructor
    {
        coot::atom_overlaps_dots_container_t aodc;
        EXPECT_EQ(aodc.dots.size(), 6u);
        EXPECT_EQ(aodc.clashes.size(), 0u);
        EXPECT_NEAR(aodc.score(), 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
