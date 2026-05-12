#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

using coot::hb_t;
using coot::HB_HYDROGEN;
using coot::HB_ACCEPTOR;
using coot::HB_DONOR;
using coot::HB_BOTH;

TEST(OracleTest, HbondTypeInfo) {
    // case: H as donor (HB_HYDROGEN) to acceptor (HB_ACCEPTOR) - H is first atom
    {
        coot::atom_overlaps_container_t::h_bond_info_t hbi(HB_HYDROGEN, HB_ACCEPTOR, "ALA", "ALA");
        
        EXPECT_TRUE(hbi.is_h_bond_H_and_acceptor);
        EXPECT_TRUE(hbi.H_is_first_atom_flag);
        EXPECT_FALSE(hbi.H_is_second_atom_flag);
        EXPECT_FALSE(hbi.is_h_bond_donor_and_acceptor);
    }
    
    // case: acceptor (HB_ACCEPTOR) to H as donor (HB_HYDROGEN) - H is second atom
    {
        coot::atom_overlaps_container_t::h_bond_info_t hbi(HB_ACCEPTOR, HB_HYDROGEN, "ALA", "ALA");
        
        EXPECT_TRUE(hbi.is_h_bond_H_and_acceptor);
        EXPECT_FALSE(hbi.H_is_first_atom_flag);
        EXPECT_TRUE(hbi.H_is_second_atom_flag);
        EXPECT_FALSE(hbi.is_h_bond_donor_and_acceptor);
    }
    
    // case: donor (HB_DONOR) to acceptor (HB_ACCEPTOR) - H-bond donor/acceptor
    {
        coot::atom_overlaps_container_t::h_bond_info_t hbi(HB_DONOR, HB_ACCEPTOR, "ALA", "ALA");
        
        EXPECT_FALSE(hbi.is_h_bond_H_and_acceptor);
        EXPECT_FALSE(hbi.H_is_first_atom_flag);
        EXPECT_FALSE(hbi.H_is_second_atom_flag);
        EXPECT_TRUE(hbi.is_h_bond_donor_and_acceptor);
    }
    
    // case: both HB_BOTH - should still be donor/acceptor
    {
        coot::atom_overlaps_container_t::h_bond_info_t hbi(HB_BOTH, HB_BOTH, "ALA", "ALA");
        
        EXPECT_FALSE(hbi.is_h_bond_H_and_acceptor);
        EXPECT_FALSE(hbi.H_is_first_atom_flag);
        EXPECT_FALSE(hbi.H_is_second_atom_flag);
        EXPECT_TRUE(hbi.is_h_bond_donor_and_acceptor);
    }
    
    // case: both HB_HYDROGEN - should NOT form H-bond
    {
        coot::atom_overlaps_container_t::h_bond_info_t hbi(HB_HYDROGEN, HB_HYDROGEN, "ALA", "ALA");
        
        EXPECT_FALSE(hbi.is_h_bond_H_and_acceptor);
        EXPECT_FALSE(hbi.H_is_first_atom_flag);
        EXPECT_FALSE(hbi.H_is_second_atom_flag);
        EXPECT_FALSE(hbi.is_h_bond_donor_and_acceptor);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}