#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_neighbor_residues_for_refinement_help) {
    // Test 1: Standard PDB
    {
        const char* pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);
        
        coot::molecule_t test_mol("test_frag", 99);
        EXPECT_EQ(test_mol.neighbouring_residues.size(), 0);
        
        test_mol.add_neighbor_residues_for_refinement_help_gemmi(st);
        
        EXPECT_EQ(test_mol.neighbouring_residues.size(), 0);
    }

    // Test 2: Empty structure (should yield 0 neighbours)
    {
        gemmi::Structure st;
        coot::molecule_t test_mol2("empty_frag", 100);
        EXPECT_EQ(test_mol2.neighbouring_residues.size(), 0);
        
        test_mol2.add_neighbor_residues_for_refinement_help_gemmi(st);
        
        EXPECT_EQ(test_mol2.neighbouring_residues.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}