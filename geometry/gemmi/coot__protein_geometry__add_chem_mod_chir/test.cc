#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>
#include "function.hh"

TEST(OracleTest, add_chem_mod_chir) {
    // Pre-populate mods to 123 entries (matching coot::protein_geometry::init_standard())
    std::map<std::string, coot::chem_mod> mods;
    for (int i = 0; i < 123; ++i) {
        mods["dummy_" + std::to_string(i)] = coot::chem_mod();
    }

    // Case 1: Valid chiral restraint with positive volume sign
    {
        size_t mods_before = mods.size();

        ChemModChirLoop loop;
        loop.tags = {"mod_id", "function", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "new_volume_sign"};
        loop.rows.push_back({"MYMOD1", "chir", "CA", "C", "N", "CB", "positive"});

        EXPECT_EQ(loop.GetLoopLength(), 1);
        EXPECT_EQ(mods_before, 123u);

        coot::protein_geometry::add_chem_mod_chir_gemmi(loop, mods);

        EXPECT_EQ(mods.size(), 124u);

        bool found = mods.count("MYMOD1") > 0;
        EXPECT_TRUE(found);

        if (found) {
            const coot::chem_mod &mod = mods.at("MYMOD1");
            EXPECT_EQ(mod.chir_mods.size(), 1u);
        }
    }

    // Case 2: Negative volume sign + multiple rows
    {
        size_t mods_before = mods.size();

        ChemModChirLoop loop2;
        loop2.tags = {"mod_id", "function", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "new_volume_sign"};
        loop2.rows.push_back({"MYMOD2", "chir", "CA", "C", "N", "CB", "negative"});
        loop2.rows.push_back({"MYMOD2", "chir", "CB", "CA", "CG", "SG", "positive"});

        EXPECT_EQ(loop2.GetLoopLength(), 2u);
        EXPECT_EQ(mods_before, 124u);

        coot::protein_geometry::add_chem_mod_chir_gemmi(loop2, mods);

        EXPECT_EQ(mods.size(), 125u);

        bool found = mods.count("MYMOD2") > 0;
        EXPECT_TRUE(found);

        if (found) {
            const coot::chem_mod &mod = mods.at("MYMOD2");
            EXPECT_EQ(mod.chir_mods.size(), 2u);
        }
    }

    // Case 3: Empty loop (edge case — should do nothing)
    {
        size_t mods_before = mods.size();

        ChemModChirLoop loop3;
        loop3.tags = {"mod_id", "function", "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3", "new_volume_sign"};

        EXPECT_EQ(loop3.GetLoopLength(), 0u);
        EXPECT_EQ(mods_before, 125u);

        coot::protein_geometry::add_chem_mod_chir_gemmi(loop3, mods);

        EXPECT_EQ(mods.size(), 125u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}