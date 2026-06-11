#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, dictionary_residue_restraints_t) {
    // Case 1: valid molecule with residues (oracle case)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        
        ASSERT_FALSE(st.models.empty());
        
        coot::dictionary_residue_restraints_t dict = 
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(&st);
        
        EXPECT_EQ(dict.bond_restraint.size(), 7u);
        EXPECT_EQ(dict.angle_restraint.size(), 7u);
        EXPECT_EQ(dict.atom_info.size(), 8u);
        EXPECT_EQ(dict.residue_info.comp_id, "MET");
        EXPECT_FALSE(dict.filled_with_bond_order_data_only_flag);
        
        EXPECT_FALSE(dict.bond_restraint.empty());
        EXPECT_EQ(dict.bond_restraint[0].atom_id_1(), " CA ");
        EXPECT_EQ(dict.bond_restraint[0].atom_id_2(), " N  ");
        
        EXPECT_FALSE(dict.angle_restraint.empty());
        EXPECT_EQ(dict.angle_restraint[0].atom_id_1(), " N  ");
        EXPECT_EQ(dict.angle_restraint[0].atom_id_2(), " CA ");
        EXPECT_EQ(dict.angle_restraint[0].atom_id_3(), " C  ");
        EXPECT_NEAR(dict.angle_restraint[0].angle(), 109.545, 1e-3);
    }
    
    // Case 2: empty structure — constructor should return with filled_with_bond_order_data_only_flag == true
    {
        gemmi::Structure empty_st;
        // empty_st has no models by default
        
        coot::dictionary_residue_restraints_t dict = 
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(&empty_st);
        
        // When no residue is found, init() is never called, so flag stays true
        EXPECT_TRUE(dict.filled_with_bond_order_data_only_flag);
        EXPECT_TRUE(dict.bond_restraint.empty());
        EXPECT_TRUE(dict.angle_restraint.empty());
        EXPECT_TRUE(dict.atom_info.empty());
    }

    // Case 3: nullptr structure
    {
        coot::dictionary_residue_restraints_t dict = 
            coot::dictionary_residue_restraints_t::dictionary_residue_restraints_t_gemmi(nullptr);
        
        EXPECT_TRUE(dict.filled_with_bond_order_data_only_flag);
        EXPECT_TRUE(dict.bond_restraint.empty());
        EXPECT_TRUE(dict.angle_restraint.empty());
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
