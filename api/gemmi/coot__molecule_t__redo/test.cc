#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, redo) {
    // Read the test PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Create molecule_t and set up its mol pointer
    coot::molecule_t mol_obj;
    mol_obj.mol = &st;
    
    // Verify initial atom count matches oracle (atom_sel starts empty)
    EXPECT_EQ(mol_obj.atom_sel.n_selected_atoms, 0);

    // ── Case 1: redo with no history (edge case — nothing to redo) ──
    {
        // Underlying redo should return nullptr (nothing to redo)
        const gemmi::Structure* mol_new = mol_obj.modification_info.redo();
        EXPECT_EQ(mol_new, nullptr);
    }

    // ── Setup: create backups so undo/redo chain exists ──
    mol_obj.modification_info.make_backup(*mol_obj.mol, "backup_A");
    mol_obj.modification_info.make_backup(*mol_obj.mol, "backup_B");

    // Verify backup count
    EXPECT_EQ(static_cast<int>(mol_obj.modification_info.save_info.size()), 2);

    // ── Undo: goes back one slot ──
    const gemmi::Structure* mol_undo = mol_obj.modification_info.undo(*mol_obj.mol);
    EXPECT_NE(mol_undo, nullptr);

    // After undo, modification_index should be 1 (oracle shows "changing to index 1")
    EXPECT_EQ(mol_obj.modification_info.modification_index, 1);

    // ── Case 2: redo with valid history ──
    {
        // Underlying redo should return non-null (valid redo available)
        const gemmi::Structure* mol_new = mol_obj.modification_info.redo();
        EXPECT_NE(mol_new, nullptr);

        // After successful redo, modification_index should advance to 2
        EXPECT_EQ(mol_obj.modification_info.modification_index, 2);
    }

    // ── Case 3: redo at end of history (nothing left to redo) ──
    {
        // Underlying redo should return nullptr (already at end)
        const gemmi::Structure* mol_new = mol_obj.modification_info.redo();
        EXPECT_EQ(mol_new, nullptr);

        // modification_index stays at 2 (no change since redo failed)
        EXPECT_EQ(mol_obj.modification_info.modification_index, 2);
    }

    // ── Verify the public redo() wrapper always returns 0 ──
    {
        int status = coot::redo_gemmi(mol_obj);
        EXPECT_EQ(status, 0);
    }

    // ── Coverage improvement: test that redo updates atom_sel on success ──
    {
        // Reset to allow a successful redo
        mol_obj.modification_info.modification_index = 0;
        // Clear atom_sel first
        mol_obj.atom_sel.clear_up();

        // Now redo should succeed and populate atom_sel
        int status = coot::redo_gemmi(mol_obj);
        EXPECT_EQ(status, 0);
        EXPECT_GT(mol_obj.atom_sel.n_selected_atoms, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
