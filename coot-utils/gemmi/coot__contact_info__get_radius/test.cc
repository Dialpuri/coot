#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, ContactInfoGetRadius) {
    coot::contact_info ci;
    
    // Populate the atom_radii vector
    ci.setup_atom_radii();
    
    // Case 1: known element (carbon)
    {
        std::string ele_C = " C";
        double r_C = ci.get_radius_gemmi(ele_C);
        EXPECT_NEAR(r_C, 0.77, 1e-4);
    }
    
    // Case 2: known element (nitrogen)
    {
        std::string ele_N = " N";
        double r_N = ci.get_radius_gemmi(ele_N);
        EXPECT_NEAR(r_N, 0.65, 1e-4);
    }
    
    // Case 3: unknown element (should return default 0.9)
    {
        std::string ele_unknown = "XY";
        double r_unknown = ci.get_radius_gemmi(ele_unknown);
        EXPECT_NEAR(r_unknown, 0.9, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}