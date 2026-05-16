#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test for coot::protein_geometry::add_chem_mod_atom_gemmi

TEST(OracleTest, add_chem_mod_atom_gemmi) {
    // case: valid _chem_mod_atom row — adds a mod atom entry
    {
        std::map<std::string, coot::chem_mod> mods;

        coot::ChemModAtomLoop loop;
        loop.tags = {
            "mod_id", "function", "atom_id",
            "new_atom_id", "new_type_symbol", "new_type_energy",
            "new_partial_charge"
        };

        // Row 0: mod_id="TESTMOD", function="add", atom_id="P",
        //         new_atom_id="PDUM", new_type_symbol="P",
        //         new_type_energy="P", new_partial_charge=0.0
        std::vector<std::string> row = {
            "TESTMOD", "add", "P",
            "PDUM", "P", "P", "0.0"
        };
        loop.rows.push_back(row);

        coot::protein_geometry::add_chem_mod_atom_gemmi(loop, mods);

        size_t mods_count = mods.size();
        EXPECT_EQ(mods_count, 1u);

        auto it = mods.find("TESTMOD");
        EXPECT_NE(it, mods.end());

        // Verify the chem_mod_atom was added correctly
        EXPECT_EQ(it->second.mod_atoms.size(), 1u);
        EXPECT_EQ(it->second.mod_atoms[0].function, "add");
        EXPECT_EQ(it->second.mod_atoms[0].atom_id, "P");
        EXPECT_EQ(it->second.mod_atoms[0].new_atom_id, "PDUM");
        EXPECT_EQ(it->second.mod_atoms[0].new_type_symbol, "P");
        EXPECT_EQ(it->second.mod_atoms[0].new_type_energy, "P");
        EXPECT_EQ(it->second.mod_atoms[0].new_partial_charge, 0.0);
    }

    // case: invalid_missing_tags — should be rejected by guard
    {
        std::map<std::string, coot::chem_mod> mods;
        coot::chem_mod tm;
        // Pre-populate to verify size doesn't change on invalid row
        tm.mod_atoms.push_back(coot::chem_mod_atom("test", "x", "y", "A", "B", 1.0));
        mods["PREEXIST"] = tm;

        coot::ChemModAtomLoop loop;
        loop.tags = {"mod_id"};

        // Only mod_id provided — other tags missing, should fail
        std::vector<std::string> row = {"BADMOD"};
        loop.rows.push_back(row);

        coot::protein_geometry::add_chem_mod_atom_gemmi(loop, mods);

        size_t mods_count = mods.size();
        EXPECT_EQ(mods_count, 1u); // should not have added BADMOD

        auto it = mods.find("BADMOD");
        EXPECT_EQ(it, mods.end());
    }

    // case: atom_id_mmdb_expand verification
    {
        EXPECT_EQ(coot::atom_id_mmdb_expand("C"), " C  ");
        EXPECT_EQ(coot::atom_id_mmdb_expand("CA"), " CA ");
        EXPECT_EQ(coot::atom_id_mmdb_expand("FE"), "FE  ");
        EXPECT_EQ(coot::atom_id_mmdb_expand("C1"), " C1 ");
        EXPECT_EQ(coot::atom_id_mmdb_expand("CCCC"), "CCCC");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}