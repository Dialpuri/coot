#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test suite for set_occupancy function
TEST(OracleTest, SetOccupancy) {
    // Setup: load structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    {  // Test case 1: Valid CID //A/10, occupancy 0.25
        std::string cid = "//A/10";
        float occ_new = 0.25f;
        
        // Get occupancy before
        float occupancy_before = -1.0f;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "10") {
                            if (!residue.atoms.empty()) {
                                occupancy_before = residue.atoms[0].occ;
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
        
        coot::molecule_t::set_occupancy_gemmi(st, cid, occ_new);
        
        // Verify occupancy was changed
        bool found = false;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "10") {
                            if (!residue.atoms.empty()) {
                                float occupancy_after = residue.atoms[0].occ;
                                EXPECT_FLOAT_EQ(occupancy_after, occ_new);
                                EXPECT_FLOAT_EQ(occupancy_before, 1.0f);
                            }
                            found = true;
                            break;
                        }
                    }
                    break;
                }
            }
        }
        if (!found) {
            FAIL() << "No atoms matched CID //A/10";
        }
    }
    
    {  // Test case 2: Another residue //A/20, occupancy 0.75
        std::string cid = "//A/20";
        float occ_new = 0.75f;
        
        // Get occupancy before
        float occupancy_before = -1.0f;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "20") {
                            if (!residue.atoms.empty()) {
                                occupancy_before = residue.atoms[0].occ;
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
        
        coot::molecule_t::set_occupancy_gemmi(st, cid, occ_new);
        
        // Verify occupancy was changed
        bool found = false;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "A") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "20") {
                            if (!residue.atoms.empty()) {
                                float occupancy_after = residue.atoms[0].occ;
                                EXPECT_FLOAT_EQ(occupancy_after, occ_new);
                            }
                            found = true;
                            break;
                        }
                    }
                    break;
                }
            }
        }
        if (!found) {
            FAIL() << "No atoms matched CID //A/20";
        }
    }
    
    {  // Test case 3: Invalid CID //Z/999, should not crash
        std::string cid = "//Z/999";
        float occ_new = 0.50f;
        
        // Should not match any atoms
        bool found = false;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "Z") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "999") {
                            found = true;
                        }
                    }
                }
            }
        }
        EXPECT_FALSE(found);
        
        // Should not crash
        coot::molecule_t::set_occupancy_gemmi(st, cid, occ_new);
        
        // Still should not match
        found = false;
        if (!st.models.empty()) {
            gemmi::Model& model = st.models[0];
            for (gemmi::Chain& chain : model.chains) {
                if (chain.name == "Z") {
                    for (gemmi::Residue& residue : chain.residues) {
                        if (std::to_string(residue.seqid.num.value) == "999") {
                            found = true;
                        }
                    }
                }
            }
        }
        EXPECT_FALSE(found);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}