#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_nbc_dist) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: standard backbone atoms, not in same residue
    {
        std::string et1 = "CA";
        std::string et2 = "N";
        bool same_res = false;
        bool same_ring = false;

        std::pair<bool, double> result = coot::get_nbc_dist_gemmi(geom, et1, et2, same_res, same_ring);

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 2.69, 1e-4);
    }

    // Case 2: same residue (gets 0.84 multiplier)
    {
        std::string et1 = "CA";
        std::string et2 = "N";
        bool same_res = true;
        bool same_ring = false;

        std::pair<bool, double> result = coot::get_nbc_dist_gemmi(geom, et1, et2, same_res, same_ring);

        EXPECT_TRUE(result.first);
        EXPECT_NEAR(result.second, 2.2596, 1e-4);
    }

    // Case 3: invalid energy types (not in energy_lib)
    {
        std::string et1 = "INVALID_TYPE_XYZ";
        std::string et2 = "ANOTHER_INVALID";
        bool same_res = false;
        bool same_ring = false;

        std::pair<bool, double> result = coot::get_nbc_dist_gemmi(geom, et1, et2, same_res, same_ring);

        EXPECT_FALSE(result.first);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}