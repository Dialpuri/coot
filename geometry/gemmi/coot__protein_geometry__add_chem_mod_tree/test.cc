#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_mod_tree) {
    coot::protein_geometry geom;
    geom.init_standard();

    // === Case 1: Valid loop with 2 rows ===
    {
        std::size_t mods_size_before = geom.mods.size();
        EXPECT_EQ(mods_size_before, 123u);

        ChemModTreeLoop loop;
        loop.tags = {
            "mod_id", "function", "atom_id", "atom_back",
            "back_type", "atom_forward", "connect_type"
        };

        // Row 0: first tree entry
        loop.rows.push_back({"MOD_A", "main", "N", "?", "?", "CA", "chain"});

        // Row 1: second tree entry for same mod
        loop.rows.push_back({"MOD_A", "branch", "CA", "N", "chain", "C", "chain"});

        int loop_len = loop.GetLoopLength();
        EXPECT_EQ(loop_len, 2);

        // Call the function
        coot::add_chem_mod_tree_gemmi(loop, geom);

        std::size_t mods_size_after = geom.mods.size();
        EXPECT_EQ(mods_size_after, 124u);
        EXPECT_TRUE(geom.mods.count("MOD_A") > 0);
        if (geom.mods.count("MOD_A")) {
            EXPECT_EQ(geom.mods["MOD_A"].tree_mods.size(), 2u);
        }
    }

    // === Case 2: Empty loop (no rows) — should add nothing ===
    {
        ChemModTreeLoop loop2;
        loop2.tags = {
            "mod_id", "function", "atom_id", "atom_back",
            "back_type", "atom_forward", "connect_type"
        };
        // No rows added

        int loop_len2 = loop2.GetLoopLength();
        EXPECT_EQ(loop_len2, 0);

        std::size_t mods_before_empty = geom.mods.size();
        coot::add_chem_mod_tree_gemmi(loop2, geom);
        std::size_t mods_after_empty = geom.mods.size();
        EXPECT_TRUE(mods_before_empty == mods_after_empty);
    }

    // === Case 3: Another mod with 1 row ===
    {
        ChemModTreeLoop loop3;
        loop3.tags = {
            "mod_id", "function", "atom_id", "atom_back",
            "back_type", "atom_forward", "connect_type"
        };

        // Valid row
        loop3.rows.push_back({"MOD_B", "main", "N", "?", "?", "CA", "chain"});

        std::size_t mods_before_modb = geom.mods.size();
        EXPECT_EQ(mods_before_modb, 124u);

        coot::add_chem_mod_tree_gemmi(loop3, geom);
        std::size_t mods_after_modb = geom.mods.size();
        EXPECT_EQ(mods_after_modb, 125u);
        EXPECT_TRUE(geom.mods.count("MOD_B") > 0);
        if (geom.mods.count("MOD_B")) {
            EXPECT_EQ(geom.mods["MOD_B"].tree_mods.size(), 1u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}