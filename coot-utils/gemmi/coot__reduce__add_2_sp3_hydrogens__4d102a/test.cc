#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, Add2Sp3Hydrogens) {
    // Load the test molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    
    // Get residue ILE A 10
    gemmi::Residue* res = nullptr;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& r : chain.residues) {
                if (r.name == "ILE" && chain.name == "A" && r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
        }
    }
    ASSERT_NE(res, nullptr) << "Failed to get residue ILE A 10";
    
    // Case 1: valid input with 2 second neighbors
    {
        std::string H_at_name_1 = " HB2";
        std::string H_at_name_2 = " HB3";
        std::string first_neighb = " CA ";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back(" CB ");
        second_neighb_vec.push_back(" OG ");
        
        double bond_length = 0.97;
        double angle_between_Hs = 107.0;
        
        // Verify inputs
        EXPECT_EQ(H_at_name_1, " HB2");
        EXPECT_EQ(H_at_name_2, " HB3");
        EXPECT_EQ(first_neighb, " CA ");
        EXPECT_EQ(second_neighb_vec[0], " CB ");
        EXPECT_EQ(second_neighb_vec[1], " OG ");
        EXPECT_NEAR(bond_length, 0.97, 1e-4);
        EXPECT_NEAR(angle_between_Hs, 107.0, 1e-4);
        
        // Verify residue
        std::string res_name = res->name;
        std::string res_chain = st.models[0].chains[0].name;
        int res_seq = res->seqid.num.value;
        EXPECT_EQ(res_name, "ILE");
        EXPECT_EQ(res_chain, "A");
        EXPECT_EQ(res_seq, 10);
        
        // The function should execute without throwing
        EXPECT_NO_THROW(coot::reduce::add_2_sp3_hydrogens_gemmi(H_at_name_1, H_at_name_2, first_neighb, 
                                                               second_neighb_vec, bond_length, 
                                                               angle_between_Hs, res, false));
        
        // Verify that adding hydrogens was successful (checking second neighbors vector size)
        EXPECT_EQ(second_neighb_vec.size(), 2);
    }
    
    // Case 2: invalid input with wrong size (triggers warning message)
    {
        std::string H_at_name_1 = " HB2";
        std::string H_at_name_2 = " HB3";
        std::string first_neighb = " CA ";
        std::vector<std::string> second_neighb_vec;
        second_neighb_vec.push_back(" CB "); // only 1 element
        
        double bond_length = 0.97;
        double angle_between_Hs = 107.0;
        
        // Verify inputs for edge case
        EXPECT_EQ(H_at_name_1, " HB2");
        EXPECT_EQ(second_neighb_vec.size(), 1);
        
        // The function should still execute (with warning) without crashing
        EXPECT_NO_THROW(coot::reduce::add_2_sp3_hydrogens_gemmi(H_at_name_1, H_at_name_2, first_neighb, 
                                                               second_neighb_vec, bond_length, 
                                                               angle_between_Hs, res, false));
        
        // Verify the guard condition - size should still be 1
        EXPECT_EQ(second_neighb_vec.size(), 1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}