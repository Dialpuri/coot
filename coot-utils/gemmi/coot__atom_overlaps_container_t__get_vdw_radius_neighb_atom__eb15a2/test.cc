#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_overlaps_container_t_get_vdw_radius_neighb_atom) {
    // Construct a minimal atom_overlaps_container_t using default constructor
    coot::atom_overlaps_container_t cont;

    // Directly populate the neighb_atom_radius vector
    cont.neighb_atom_radius.clear();
    cont.neighb_atom_radius.push_back(1.70);  // idx 0: typical C vdw radius
    cont.neighb_atom_radius.push_back(1.50);  // idx 1: typical N vdw radius
    cont.neighb_atom_radius.push_back(1.40);  // idx 2: typical O vdw radius

    // Case 1: valid index 0
    {
        int idx = 0;
        double r = cont.get_vdw_radius_neighb_atom_gemmi(idx);
        EXPECT_NEAR(r, 1.7, 1e-4);
    }

    // Case 2: valid index 1
    {
        int idx = 1;
        double r = cont.get_vdw_radius_neighb_atom_gemmi(idx);
        EXPECT_NEAR(r, 1.5, 1e-4);
    }

    // Case 3: last valid index 2
    {
        int idx = 2;
        double r = cont.get_vdw_radius_neighb_atom_gemmi(idx);
        EXPECT_NEAR(r, 1.4, 1e-4);
    }

    // Case 4: out-of-bounds index (no index checking in the function)
    {
        int idx = 5;  // beyond the vector
        EXPECT_NO_THROW(cont.get_vdw_radius_neighb_atom_gemmi(idx));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
