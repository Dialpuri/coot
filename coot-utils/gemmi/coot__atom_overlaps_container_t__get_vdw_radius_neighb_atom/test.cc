#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, AtomOverlapsVDWRadius) {
    coot::atom_overlaps_container_t aoc;
    aoc.neighb_atom_radius.push_back(1.55);
    aoc.neighb_atom_radius.push_back(1.52);
    
    {
        double r0 = aoc.get_vdw_radius_neighb_atom_gemmi(0);
        EXPECT_NEAR(r0, 1.55, 1e-4);
    }
    
    {
        int last_idx = aoc.neighb_atom_radius.size() - 1;
        double r_last = aoc.get_vdw_radius_neighb_atom_gemmi(last_idx);
        EXPECT_NEAR(r_last, 1.52, 1e-4);
    }
    
    {
        int bad_idx = aoc.neighb_atom_radius.size() + 100;
        double r_bad = aoc.get_vdw_radius_neighb_atom_gemmi(bad_idx);
        EXPECT_NEAR(r_bad, 1.52, 1e-4);
    }
    
    {
        int neg_idx = -1;
        double r_neg = aoc.get_vdw_radius_neighb_atom_gemmi(neg_idx);
        EXPECT_NEAR(r_neg, 0.0, 1e-10);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}