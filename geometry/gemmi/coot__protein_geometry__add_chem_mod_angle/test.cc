#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_mod_angle) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Valid chem_mod_angle loop with one row
    {
        std::vector<std::string> tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2", "atom_id_3",
            "new_value_angle", "new_value_angle_esd"
        };

        std::vector<std::vector<std::string>> rows;
        rows.push_back({"TEST_MOD", "angle", "N", "CA", "C", "120.0", "5.0"});

        std::size_t mods_before = geom.mods.size();
        EXPECT_EQ(mods_before, 123u);

        coot::add_chem_mod_angle_gemmi(tags, rows, geom.mods);

        EXPECT_EQ(geom.mods.size(), 124u);

        auto it = geom.mods.find("TEST_MOD");
        EXPECT_NE(it, geom.mods.end());
        EXPECT_EQ(it->second.angle_mods.size(), 1u);

        const coot::chem_mod_angle &cma = it->second.angle_mods[0];
        EXPECT_EQ(static_cast<int>(cma.function), 0);
        EXPECT_EQ(cma.atom_id_1, " N  ");
        EXPECT_EQ(cma.atom_id_2, " CA ");
        EXPECT_EQ(cma.atom_id_3, " C  ");
        EXPECT_NEAR(cma.new_value_angle, 120.0, 1e-4);
        EXPECT_NEAR(cma.new_value_angle_esd, 5.0, 1e-4);
    }

    // Case 2: Empty loop — should add nothing
    {
        std::vector<std::string> tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2", "atom_id_3",
            "new_value_angle", "new_value_angle_esd"
        };

        std::vector<std::vector<std::string>> rows_empty;

        std::size_t mods_before = geom.mods.size();
        EXPECT_EQ(mods_before, 124u);

        coot::add_chem_mod_angle_gemmi(tags, rows_empty, geom.mods);

        EXPECT_EQ(geom.mods.size(), 124u);
    }

    // Case 3: Loop with two rows for same mod_id
    {
        std::vector<std::string> tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2", "atom_id_3",
            "new_value_angle", "new_value_angle_esd"
        };

        std::vector<std::vector<std::string>> rows;

        // Row 0
        rows.push_back({"TEST_MOD2", "angle", "O", "CB", "CG", "110.5", "2.0"});

        // Row 1
        rows.push_back({"TEST_MOD2", "angle", "N", "CA", "CB", "109.0", "1.5"});

        std::size_t mods_before = geom.mods.size();
        EXPECT_EQ(mods_before, 124u);

        coot::add_chem_mod_angle_gemmi(tags, rows, geom.mods);

        EXPECT_EQ(geom.mods.size(), 125u);

        auto it = geom.mods.find("TEST_MOD2");
        EXPECT_NE(it, geom.mods.end());
        EXPECT_EQ(it->second.angle_mods.size(), 2u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
