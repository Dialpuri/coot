#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, MatchContainerAdd) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Get atoms from residue 10 (same residue)
    gemmi::Atom* at1 = nullptr;
    gemmi::Atom* at2 = nullptr;
    
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == 10) {
                for (gemmi::Atom& atom : res.atoms) {
                    if (atom.name == "CA" && !at1) at1 = &atom;
                    if (atom.name == "C" && !at2) at2 = &atom;
                }
            }
        }
    }
    
    ASSERT_NE(at1, nullptr) << "Failed to get atom at1 (CA)";
    ASSERT_NE(at2, nullptr) << "Failed to get atom at2 (C)";
    
    // Test case 1: add two atoms that belong to the same residue
    {
        coot::match_container_t match_set;
        EXPECT_EQ(match_set.matches.size(), 0u);
        
        match_set.add_gemmi(at1, at2, model);
        EXPECT_EQ(match_set.matches.size(), 1u);
        
        if (match_set.matches.size() > 0) {
            EXPECT_EQ(match_set.matches[0].atom_pairs.size(), 1u);
        }
    }
    
    // Test case 2: add atoms from different residue pairs
    {
        gemmi::Atom* at3 = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == 11) {
                    for (gemmi::Atom& atom : res.atoms) {
                        if (atom.name == "CA" && !at3) at3 = &atom;
                    }
                }
            }
        }
        
        if (at3) {
            coot::match_container_t match_set;
            match_set.add_gemmi(at1, at3, model); // at1 is from residue 10, at3 is from residue 11
            
            EXPECT_EQ(match_set.matches.size(), 1u);
            if (match_set.matches.size() > 0) {
                EXPECT_EQ(match_set.matches[0].atom_pairs.size(), 1u);
            }
        }
    }
    
    // Test case 3: add same atom pair twice (should group, not duplicate)
    {
        coot::match_container_t match_set;
        match_set.add_gemmi(at1, at2, model);
        match_set.add_gemmi(at1, at2, model); // add again to same residue pair
        
        EXPECT_EQ(match_set.matches.size(), 1u);
        if (match_set.matches.size() > 0) {
            EXPECT_EQ(match_set.matches[0].atom_pairs.size(), 2u);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}