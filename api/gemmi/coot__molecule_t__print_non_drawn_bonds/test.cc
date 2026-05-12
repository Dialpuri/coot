#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, PrintNonDrawnBonds) {
    // Create molecule from PDB file
    coot::molecule_t mol("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test 1: Initial state (no bonds set)
    EXPECT_EQ(mol.no_bonds_to_these_atom_indices.size(), 0);
    mol.print_non_drawn_bonds();
    
    // Test 2: Add first atom index (index 0)
    mol.no_bonds_to_these_atom_indices.insert(0);
    EXPECT_EQ(mol.no_bonds_to_these_atom_indices.size(), 1);
    mol.print_non_drawn_bonds();
    
    // Test 3: Add multiple atom indices
    mol.no_bonds_to_these_atom_indices.insert(5);
    mol.no_bonds_to_these_atom_indices.insert(10);
    EXPECT_EQ(mol.no_bonds_to_these_atom_indices.size(), 3);
    mol.print_non_drawn_bonds();
    
    // Test 4: Add an invalid atom index (out of range)
    mol.no_bonds_to_these_atom_indices.insert(9999);
    EXPECT_EQ(mol.no_bonds_to_these_atom_indices.size(), 4);
    mol.print_non_drawn_bonds();
    
    // Test 5: Clear and test empty set
    mol.no_bonds_to_these_atom_indices.clear();
    EXPECT_EQ(mol.no_bonds_to_these_atom_indices.size(), 0);
    mol.print_non_drawn_bonds();
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}