#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, add_chem_mod_plane) {
  // ================================================================
  // Case 1: Oracle's original case — loop with 3 rows for ILE/PL1
  // Expected: function does nothing (ierr_tot == 5 per iteration)
  // ================================================================
  {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::cif::Loop loop;
    loop.tags = {
        "_chem_mod_plane_atom.mod_id",
        "_chem_mod_plane_atom.plane_id",
        "_chem_mod_plane_atom.function",
        "_chem_mod_plane_atom.atom_id",
        "_chem_mod_plane_atom.new_dist_esd"
    };

    loop.add_row({"ILE", "PL1", "add", "CA", "0.01"});
    loop.add_row({"ILE", "PL1", "add", "C", "0.02"});
    loop.add_row({"ILE", "PL1", "add", "CB", "0.015"});

    EXPECT_EQ(loop.length(), 3);

    size_t mods_size_before = geom.mods.size();
    EXPECT_EQ(mods_size_before, 123u);

    // Call the function
    coot::add_chem_mod_plane_gemmi(geom, loop);

    // After: size unchanged, ILE not found in mods
    EXPECT_EQ(geom.mods.size(), 123u);
    EXPECT_FALSE(geom.mods.count("ILE"));
  }

  // ================================================================
  // Case 2: Complementary — empty loop (0 rows) should be handled
  // without crash / no mutation
  // ================================================================
  {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::cif::Loop loop;
    // tags is empty, so loop.length() would divide by zero — check tags instead
    EXPECT_TRUE(loop.tags.empty());
    size_t mods_size_before = geom.mods.size();

    EXPECT_NO_THROW(coot::add_chem_mod_plane_gemmi(geom, loop));

    EXPECT_EQ(geom.mods.size(), mods_size_before);
  }

  // ================================================================
  // Case 3: Complementary — 1-row loop with "remove" function for
  // a mod that doesn't exist yet (edge: single iteration, remove)
  // ================================================================
  {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::cif::Loop loop;
    loop.tags = {
        "_chem_mod_plane_atom.mod_id",
        "_chem_mod_plane_atom.plane_id",
        "_chem_mod_plane_atom.function",
        "_chem_mod_plane_atom.atom_id",
        "_chem_mod_plane_atom.new_dist_esd"
    };

    // Single row: try to remove a plane for a non-existent mod
    loop.add_row({"GLY", "PL9", "remove", "N", "0.0"});

    EXPECT_EQ(loop.length(), 1);

    size_t mods_size_before = geom.mods.size();
    EXPECT_NO_THROW(coot::add_chem_mod_plane_gemmi(geom, loop));
    EXPECT_EQ(geom.mods.size(), mods_size_before);
    EXPECT_FALSE(geom.mods.count("GLY"));
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}