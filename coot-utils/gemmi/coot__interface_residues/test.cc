#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(InterfaceResiduesTest, InterfaceResidues) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Test 1: Try chains A and B - check how many residues are in each chain
    {
        int n_residues_A = 0;
        int n_residues_B = 0;
        for (gemmi::Chain& chain : st.models[0].chains) {
            if (chain.name == "A") n_residues_A = chain.residues.size();
            if (chain.name == "B") n_residues_B = chain.residues.size();
        }
        EXPECT_EQ(n_residues_A, 267) << "Chain A should have 267 residues";
        EXPECT_EQ(n_residues_B, 1) << "Chain B should have 1 residue";
        
        float min_dist = 3.5;
        std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>> result = 
            coot::interface_residues_gemmi(st, "A", "B", min_dist);
        
        EXPECT_EQ(result.first.size(), 2) << "Chain A interface residues should be 2";
        EXPECT_EQ(result.second.size(), 1) << "Chain B interface residues should be 1";
    }
    
    // Test 2: Try with a non-existent chain
    {
        std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>> result = 
            coot::interface_residues_gemmi(st, "C", "A", 3.5);
        
        EXPECT_EQ(result.first.size(), 0) << "Non-existent chain C should have 0 interface residues";
        EXPECT_EQ(result.second.size(), 2) << "Chain A interface with non-existent chain should be 2";
    }
    
    // Test 3: Try with very small distance (no contacts expected)
    {
        std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>> result = 
            coot::interface_residues_gemmi(st, "A", "A", 0.5);
        
        EXPECT_EQ(result.first.size(), 0) << "Very small distance should have 0 interface residues";
        EXPECT_EQ(result.second.size(), 0) << "Very small distance should have 0 interface residues";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}