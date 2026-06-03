#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_mod_tor) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Loop with one torsion row
    {
        gemmi::cif::Loop loop;

        // Add column tags (required before adding rows)
        loop.tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2",
            "atom_id_3", "atom_id_4", "new_value_angle",
            "new_value_angle_esd", "new_period"
        };

        // Row 0: a simple dihedral
        loop.add_row({"TEST_MOD", "angle", "C", "N", "CA", "CB", "180.0", "10.0", "4"});

        EXPECT_EQ(geom.mods.size(), 123u);

        add_chem_mod_tor_gemmi(geom, loop);

        EXPECT_EQ(geom.mods.size(), 124u);

        auto it = geom.mods.find("TEST_MOD");
        EXPECT_NE(it, geom.mods.end());
        {
            coot::chem_mod &cm = it->second;
            EXPECT_EQ(it->first, "TEST_MOD");
            EXPECT_EQ(cm.tor_mods.size(), 1u);
            {
                const coot::chem_mod_tor &tor = cm.tor_mods[0];
                EXPECT_EQ(tor.atom_id_1, " C  ");
                EXPECT_EQ(tor.atom_id_2, " N  ");
                EXPECT_EQ(tor.atom_id_3, " CA ");
                EXPECT_EQ(tor.atom_id_4, " CB ");
                EXPECT_NEAR(tor.new_value_angle, 180.0, 1e-4);
                EXPECT_NEAR(tor.new_value_angle_esd, 10.0, 1e-4);
                EXPECT_EQ(tor.new_period, 4);
            }
        }
    }

    // Case 2: Empty loop — should be a no-op
    {
        gemmi::cif::Loop empty_loop;
        empty_loop.tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2",
            "atom_id_3", "atom_id_4", "new_value_angle",
            "new_value_angle_esd", "new_period"
        };

        size_t before_size = geom.mods.size();
        EXPECT_EQ(before_size, 124u);

        add_chem_mod_tor_gemmi(geom, empty_loop);

        size_t after_size = geom.mods.size();
        EXPECT_EQ(after_size, 124u);
    }

    // Case 3: Loop with bad data (missing field) — should be skipped
    {
        gemmi::cif::Loop bad_loop;
        bad_loop.tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2",
            "atom_id_3", "atom_id_4", "new_value_angle",
            "new_value_angle_esd", "new_period"
        };

        // Only fill some fields — leave new_value_angle and new_period empty (use "." for null)
        bad_loop.add_row({"BAD_MOD", "angle", "C", "N", "CA", "CB", ".", ".", "."});

        size_t before_size = geom.mods.size();
        EXPECT_EQ(before_size, 124u);

        add_chem_mod_tor_gemmi(geom, bad_loop);

        size_t after_size = geom.mods.size();
        EXPECT_EQ(after_size, 124u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}