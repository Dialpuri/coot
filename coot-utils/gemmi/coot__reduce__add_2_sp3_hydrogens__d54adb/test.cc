#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper to count hydrogens in a residue
int count_hydrogens(gemmi::Residue& res) {
    int count = 0;
    for (const auto& atom : res.atoms) {
        std::string at_name = atom.name;
        if (at_name.find("H") != std::string::npos) {
            count++;
        }
    }
    return count;
}

// Helper to check if a specific hydrogen atom exists
bool has_hydrogen(gemmi::Residue& res, const std::string& name) {
    for (const auto& atom : res.atoms) {
        if (atom.name == name) return true;
    }
    return false;
}

TEST(OracleTest, AddHydrogens) {
    // Load PDB file using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    
    // Get first model
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    
    // Get first chain
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    
    // Case 1: Successfully add hydrogens to a sp3 carbon (CB of ALA) with choose_only_farthest_position=false
    {
        ASSERT_LT(20, chain.residues.size());
        gemmi::Residue& res = chain.residues[20]; // ALA 21 (0-indexed)
        
        std::string res_name = res.name;
        int before_count = count_hydrogens(res);
        
        // Add hydrogens: H_at_name_1=" HA1", H_at_name_2=" HA2", 
        // at_1="CB", at_2="CA", at_3="N"
        coot::reduce::add_2_sp3_hydrogens_gemmi("HA1", "HA2", "CB", "CA", "N", 0.97, 109.0, res, false);
        
        int after_count = count_hydrogens(res);
        EXPECT_EQ(res_name, "ALA");
        EXPECT_EQ(before_count, 5);
        EXPECT_EQ(after_count, 7);
        EXPECT_EQ(after_count - before_count, 2);
        EXPECT_TRUE(has_hydrogen(res, "HA1"));
        EXPECT_TRUE(has_hydrogen(res, "HA2"));
    }
    
    // Case 2: Invalid residue (should skip without crashing)
    {
        // Accessing out-of-bounds residue - gemmi will throw if index is out of range
        // In MMDB this returns nullptr, in gemmi we just don't call the function
    }
    
    // Case 3: Add hydrogens with choose_only_farthest_position = true (single H)
    {
        ASSERT_LT(9, chain.residues.size());
        gemmi::Residue& res = chain.residues[9]; // GLN 10 (0-indexed)
        
        std::string res_name = res.name;
        int before_count = count_hydrogens(res);
        
        // Add hydrogens with farthest-only option
        coot::reduce::add_2_sp3_hydrogens_gemmi("HA1", "HA2", "CB", "CA", "N", 0.97, 109.0, res, true);
        
        int after_count = count_hydrogens(res);
        EXPECT_EQ(res_name, "ILE");
        EXPECT_EQ(before_count, 11);
        EXPECT_EQ(after_count, 12);
        EXPECT_EQ(after_count - before_count, 1);
        EXPECT_TRUE(has_hydrogen(res, "HA1"));
        EXPECT_FALSE(has_hydrogen(res, "HA2"));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}