#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, get_molecule_selection_as_json_gemmi) {
    // Load the PDB file
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB";
    
    // Test 1: Valid CID - chain A, residue 10 (ILE)
    {
        std::string cid1 = "//A/10";
        std::string result1 = coot::get_molecule_selection_as_json_gemmi(st, cid1);
        
        EXPECT_EQ(result1.size(), 2453u);
        EXPECT_FALSE(result1.empty());
        EXPECT_NE(result1.find("ILE"), std::string::npos);
        EXPECT_NE(result1.find("seqnum"), std::string::npos);
        EXPECT_NE(result1.find("chain_id"), std::string::npos);
    }
    
    // Test 2: Valid CID - chain A, residue 1 (MET)
    {
        std::string cid2 = "//A/1";
        std::string result2 = coot::get_molecule_selection_as_json_gemmi(st, cid2);
        
        EXPECT_EQ(result2.size(), 2440u);
        EXPECT_FALSE(result2.empty());
        EXPECT_NE(result2.find("MET"), std::string::npos);
    }
    
    // Test 3: Valid CID - chain B, residue 1 (HOH)
    {
        std::string cid3 = "//B/1";
        std::string result3 = coot::get_molecule_selection_as_json_gemmi(st, cid3);
        
        EXPECT_EQ(result3.size(), 564u);
        EXPECT_FALSE(result3.empty());
        EXPECT_NE(result3.find("HOH"), std::string::npos);
    }
    
    // Test 4: Invalid CID - non-existent chain/residue (should return empty string)
    {
        std::string cid4 = "//Z/999";
        std::string result4 = coot::get_molecule_selection_as_json_gemmi(st, cid4);
        
        EXPECT_EQ(result4.size(), 0u);
        EXPECT_TRUE(result4.empty());
    }
    
    // Test 5: Valid CID - chain A, residue 2 (GLU)
    {
        std::string cid5 = "//A/2";
        std::string result5 = coot::get_molecule_selection_as_json_gemmi(st, cid5);
        
        EXPECT_EQ(result5.size(), 2710u);
        EXPECT_FALSE(result5.empty());
        EXPECT_NE(result5.find("GLU"), std::string::npos);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
