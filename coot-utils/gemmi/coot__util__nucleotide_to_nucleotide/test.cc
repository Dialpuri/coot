#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <iostream>
#include <string>
#include "function.hh"

// Helper function to create a residue for testing
static gemmi::Residue create_test_residue(const std::string& res_name) {
    gemmi::Residue res;
    res.name = res_name;
    
    // Create a simple set of atoms for testing
    // Adenine/Guanine purine core
    if (res_name == "A" || res_name == "G") {
        gemmi::Atom n9, c8, n7, c5, c4, n1, c2, n3, c6;
        n9.name = " N9 "; n9.pos = gemmi::Position(0, 0, 0);
        c8.name = " C8 "; c8.pos = gemmi::Position(1, 0, 0);
        n7.name = " N7 "; n7.pos = gemmi::Position(1, 1, 0);
        c5.name = " C5 "; c5.pos = gemmi::Position(0, 1, 0);
        c4.name = " C4 "; c4.pos = gemmi::Position(-1, 1, 0);
        n1.name = " N1 "; n1.pos = gemmi::Position(-1, 0, 0);
        c2.name = " C2 "; c2.pos = gemmi::Position(-1, -1, 0);
        n3.name = " N3 "; n3.pos = gemmi::Position(0, -1, 0);
        c6.name = " C6 "; c6.pos = gemmi::Position(1, -1, 0);
        
        res.atoms.push_back(n9);
        res.atoms.push_back(c8);
        res.atoms.push_back(n7);
        res.atoms.push_back(c5);
        res.atoms.push_back(c4);
        res.atoms.push_back(n1);
        res.atoms.push_back(c2);
        res.atoms.push_back(n3);
        res.atoms.push_back(c6);
        
        if (res_name == "A") {
            gemmi::Atom n6;
            n6.name = " N6 "; n6.pos = gemmi::Position(2, -1, 0);
            res.atoms.push_back(n6);
        } else if (res_name == "G") {
            gemmi::Atom o6, n2;
            o6.name = " O6 "; o6.pos = gemmi::Position(2, -1, 0);
            n2.name = " N2 "; n2.pos = gemmi::Position(-2, 0, 0);
            res.atoms.push_back(o6);
            res.atoms.push_back(n2);
        }
    } else if (res_name == "U" || res_name == "C" || res_name == "T") {
        // Pyrimidine ring
        gemmi::Atom n1, c2, n3, c4, c5, c6;
        n1.name = " N1 "; n1.pos = gemmi::Position(0, 0, 0);
        c2.name = " C2 "; c2.pos = gemmi::Position(1, 0, 0);
        n3.name = " N3 "; n3.pos = gemmi::Position(1, 1, 0);
        c4.name = " C4 "; c4.pos = gemmi::Position(0, 1, 0);
        c5.name = " C5 "; c5.pos = gemmi::Position(-1, 1, 0);
        c6.name = " C6 "; c6.pos = gemmi::Position(-1, 0, 0);
        
        res.atoms.push_back(n1);
        res.atoms.push_back(c2);
        res.atoms.push_back(n3);
        res.atoms.push_back(c4);
        res.atoms.push_back(c5);
        res.atoms.push_back(c6);
        
        if (res_name == "U") {
            gemmi::Atom o2, o4;
            o2.name = " O2 "; o2.pos = gemmi::Position(2, 0, 0);
            o4.name = " O4 "; o4.pos = gemmi::Position(0, 2, 0);
            res.atoms.push_back(o2);
            res.atoms.push_back(o4);
        } else if (res_name == "C") {
            gemmi::Atom o2, n4;
            o2.name = " O2 "; o2.pos = gemmi::Position(2, 0, 0);
            n4.name = " N4 "; n4.pos = gemmi::Position(0, 2, 0);
            res.atoms.push_back(o2);
            res.atoms.push_back(n4);
        } else if (res_name == "T") {
            gemmi::Atom o2, o4, c7;
            o2.name = " O2 "; o2.pos = gemmi::Position(2, 0, 0);
            o4.name = " O4 "; o4.pos = gemmi::Position(0, 2, 0);
            c7.name = " C7 "; c7.pos = gemmi::Position(2, 1, 0);
            res.atoms.push_back(o2);
            res.atoms.push_back(o4);
            res.atoms.push_back(c7);
        }
    }
    
    return res;
}

TEST(OracleTest, nucleotide_to_nucleotide) {
    // Test case 1: Adenine to Adenine (should succeed)
    {
        gemmi::Residue adenine = create_test_residue("A");
        bool use_old_names = false;
        std::pair<bool, clipper::RTop_orth> rtop_pair =
            coot::util::nucleotide_to_nucleotide_gemmi(adenine, adenine, use_old_names);
        
        EXPECT_TRUE(rtop_pair.first);
        if (rtop_pair.first) {
            // rtop is defined
        }
    }
    
    // Test case 2: Adenine to Guanine (different bases, should succeed with alignment)
    {
        gemmi::Residue adenine = create_test_residue("A");
        gemmi::Residue guanine = create_test_residue("G");
        bool use_old_names = false;
        std::pair<bool, clipper::RTop_orth> rtop_pair =
            coot::util::nucleotide_to_nucleotide_gemmi(adenine, guanine, use_old_names);
        
        EXPECT_TRUE(rtop_pair.first);
        if (rtop_pair.first) {
            // rtop is defined
        }
    }
    
    // Test case 3: Guanine to Adenine (reverse direction)
    {
        gemmi::Residue guanine = create_test_residue("G");
        gemmi::Residue adenine = create_test_residue("A");
        bool use_old_names = false;
        std::pair<bool, clipper::RTop_orth> rtop_pair =
            coot::util::nucleotide_to_nucleotide_gemmi(guanine, adenine, use_old_names);
        
        EXPECT_TRUE(rtop_pair.first);
        if (rtop_pair.first) {
            // rtop is defined
        }
    }
    
    // Test case 4: Uracil to Uracil (additional coverage)
    {
        gemmi::Residue uracil = create_test_residue("U");
        bool use_old_names = false;
        std::pair<bool, clipper::RTop_orth> rtop_pair =
            coot::util::nucleotide_to_nucleotide_gemmi(uracil, uracil, use_old_names);
        
        EXPECT_TRUE(rtop_pair.first);
        if (rtop_pair.first) {
            // rtop is defined
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}