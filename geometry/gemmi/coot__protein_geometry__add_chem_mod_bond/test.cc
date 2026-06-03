#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_mod_bond_gemmi) {
    // Case 1: normal loop with two rows
    {
        std::map<std::string, coot::chem_mod> mods;

        std::vector<std::string> tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2",
            "new_type", "new_value_dist", "new_value_dist_esd"
        };

        std::vector<std::vector<std::string>> rows;

        // Row 0
        rows.push_back({"MYMOD", "add", "SG", "SG", "2", "2.0", "0.02"});

        // Row 1
        rows.push_back({"MYMOD", "delete", "OG", "HG", "4", "1.8", "0.01"});

        EXPECT_EQ(rows.size(), size_t(2));

        EXPECT_EQ(tags[0], "mod_id");
        EXPECT_EQ(rows[0][0], "MYMOD");
        EXPECT_EQ(rows[1][0], "MYMOD");

        EXPECT_EQ(mods.size(), size_t(0));

        coot::protein_geometry::add_chem_mod_bond_gemmi(rows, tags, mods);

        EXPECT_EQ(mods.size(), size_t(1));

        auto it = mods.find("MYMOD");
        EXPECT_NE(it, mods.end());
        {
            coot::chem_mod mod = it->second;
            EXPECT_EQ(mod.bond_mods.size(), size_t(2));

            // Bond 0
            coot::chem_mod_bond const& cmb0 = mod.bond_mods[0];
            EXPECT_EQ(cmb0.function, int(0));  // CHEM_MOD_FUNCTION_ADD
            EXPECT_EQ(cmb0.atom_id_1, " SG ");
            EXPECT_EQ(cmb0.atom_id_2, " SG ");
            EXPECT_EQ(cmb0.new_type, "2");
            EXPECT_NEAR(cmb0.new_value_dist, 2.0, 1e-4);

            // Bond 1
            coot::chem_mod_bond const& cmb1 = mod.bond_mods[1];
            EXPECT_EQ(cmb1.function, int(1));  // CHEM_MOD_FUNCTION_DELETE
            EXPECT_EQ(cmb1.atom_id_1, " OG ");
            EXPECT_EQ(cmb1.atom_id_2, "HG  ");
            EXPECT_EQ(cmb1.new_type, "4");
            EXPECT_NEAR(cmb1.new_value_dist, 1.8, 1e-4);
        }
    }

    // Case 2: empty loop — should do nothing
    {
        std::map<std::string, coot::chem_mod> mods;

        std::vector<std::string> tags = {
            "mod_id", "function", "atom_id_1", "atom_id_2",
            "new_type", "new_value_dist", "new_value_dist_esd"
        };

        std::vector<std::vector<std::string>> rows;
        rows.push_back({"MYMOD", "add", "SG", "SG", "2", "2.0", "0.02"});
        rows.push_back({"MYMOD", "delete", "OG", "HG", "4", "1.8", "0.01"});

        coot::protein_geometry::add_chem_mod_bond_gemmi(rows, tags, mods);

        // Now test empty loop
        std::vector<std::vector<std::string>> rows_empty;
        EXPECT_EQ(rows_empty.size(), size_t(0));

        EXPECT_EQ(mods.size(), size_t(1));
        coot::protein_geometry::add_chem_mod_bond_gemmi(rows_empty, tags, mods);
        EXPECT_EQ(mods.size(), size_t(1));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}