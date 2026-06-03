#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cifdoc.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, add_chem_mod_angle) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case: normal loop with 2 angle rows
    {
        EXPECT_EQ(geom.mods.size(), size_t(123));

        gemmi::cif::Loop loop;
        loop.tags = {"mod_id", "function", "atom_id_1", "atom_id_2", "atom_id_3", "new_value_angle", "new_value_angle_esd"};

        // Row 0
        loop.add_row({"MOD1", "angle", "N", "CA", "C", "120.0", "2.0"});
        // Row 1
        loop.add_row({"MOD1", "angle", "CA", "C", "CB", "110.0", "1.5"});

        EXPECT_EQ(loop.length(), 2u);
        add_chem_mod_angle_gemmi(geom, loop);

        EXPECT_EQ(geom.mods.size(), size_t(124));

        auto it = geom.mods.find("MOD1");
        ASSERT_NE(it, geom.mods.end());
        EXPECT_EQ(it->second.angle_mods.size(), size_t(2));

        // Angle 0
        coot::chem_mod_angle const &cma0 = it->second.angle_mods[0];
        EXPECT_EQ(static_cast<int>(cma0.function), 0);
        EXPECT_EQ(cma0.atom_id_1, " N  ");
        EXPECT_EQ(cma0.atom_id_2, " CA ");
        EXPECT_EQ(cma0.atom_id_3, " C  ");
        EXPECT_NEAR(cma0.new_value_angle, 120.0, 1e-4);
        EXPECT_NEAR(cma0.new_value_angle_esd, 2.0, 1e-4);

        // Angle 1
        coot::chem_mod_angle const &cma1 = it->second.angle_mods[1];
        EXPECT_EQ(static_cast<int>(cma1.function), 0);
        EXPECT_EQ(cma1.atom_id_1, " CA ");
        EXPECT_EQ(cma1.atom_id_2, " C  ");
        EXPECT_EQ(cma1.atom_id_3, " CB ");
        EXPECT_NEAR(cma1.new_value_angle, 110.0, 1e-4);
        EXPECT_NEAR(cma1.new_value_angle_esd, 1.5, 1e-4);
    }

    // Case: empty loop (0 rows)
    {
        gemmi::cif::Loop empty_loop;
        empty_loop.tags = {"mod_id", "function", "atom_id_1", "atom_id_2", "atom_id_3", "new_value_angle", "new_value_angle_esd"};

        EXPECT_EQ(empty_loop.length(), 0u);
        EXPECT_EQ(geom.mods.size(), size_t(124));
        EXPECT_NO_THROW(add_chem_mod_angle_gemmi(geom, empty_loop));
        EXPECT_EQ(geom.mods.size(), size_t(124));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}