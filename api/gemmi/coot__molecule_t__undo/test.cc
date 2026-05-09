#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <iostream>
#include "function.hh"

// Test the undo functionality with gemmi
TEST(UndoTest, SimpleUndo) {
    // Read a PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    coot::molecule_t mol;
    mol.atom_sel = st;
    mol.really_do_backups = true;
    
    // Check initial state
    EXPECT_TRUE(mol.really_do_backups);
    EXPECT_EQ(mol.modification_info.save_info.size(), 0);
    
    // Make a backup first - this is critical for undo to work
    std::string backup_result = mol.modification_info.make_backup(st, "initial state");
    EXPECT_TRUE(backup_result.empty()) << "Backup failed: " << backup_result;
    EXPECT_EQ(mol.modification_info.save_info.size(), 1);
    
    // Now call undo
    int undo_result = mol.undo();
    EXPECT_EQ(undo_result, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}