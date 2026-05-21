#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include "function.hh"

TEST(OracleTest, torsion_general_change_by) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get model, chain A, residue 15 (TYR)
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& res = chain.residues[14];  // 0-indexed, so residue 14 is seq 15
    
    // Build list of all coordinates and find atom indices
    std::vector<gemmi::Position> coords;
    int atom_N_idx = -1, atom_CA_idx = -1, atom_CB_idx = -1, atom_CG_idx = -1;
    for (size_t i = 0; i < res.atoms.size(); ++i) {
        const std::string& name = res.atoms[i].name;
        coords.push_back(res.atoms[i].pos);
        // Check with trimmed name (gemmi strips trailing spaces)
        std::string trimmed = name;
        while (!trimmed.empty() && trimmed.back() == ' ') trimmed.pop_back();
        if (trimmed == "N") atom_N_idx = static_cast<int>(i);
        else if (trimmed == "CA") atom_CA_idx = static_cast<int>(i);
        else if (trimmed == "CB") atom_CB_idx = static_cast<int>(i);
        else if (trimmed == "CG") atom_CG_idx = static_cast<int>(i);
    }
    
    // From original code:
    // clicked_atom_indices[0] = 0 (N) 
    // clicked_atom_indices[1] = 1 (CA)
    // clicked_atom_indices[2] = 4 (CB)
    // clicked_atom_indices[3] = 5 (CG)
    // 
    // Rotation about bond between atoms[1] and atoms[4] (CA-CB)
    // Rotating subtree of atoms[4] (CB and its children: CG, CD1, CD2, CE1, CE2, CZ, OH)
    
    // Children of atom CB (index 4) are: CG (5), CD1 (6), CD2 (7), CE1 (8), CE2 (9), CZ (10), OH (11)
    std::vector<int> children_of_CB = {5, 6, 7, 8, 9, 10, 11};
    
    // Case 1: valid diff of 30 degrees
    {
        // BEFORE coordinates (original PDB coordinates, unrotated)
        // These are the coordinates from the PDB file before any rotation
        double before_coords[12][3] = {
            {17.343, 12.361, 58.028},   // N (index 0)
            {16.400, 11.464, 58.677},    // CA (index 1)
            {16.184, 11.852, 60.139},  // C (index 2)
            {15.067, 12.195, 60.529},  // O (index 3)
            {16.830, 10.003, 58.517},   // CB (index 4)
            {16.483, 9.391, 57.182},   // CG (index 5) - ORIGINAL unrotated
            {15.177, 9.399, 56.709},   // CD1 (index 6)
            {17.456, 8.794, 56.394},   // CD2 (index 7)
            {14.847, 8.838, 55.485},   // CE1 (index 8)
            {17.143, 8.226, 55.168},   // CE2 (index 9)
            {15.834, 8.248, 54.713},   // CZ (index 10)
            {15.518, 7.688, 53.499},   // OH (index 11)
        };
        
        // Verify BEFORE coordinates match PDB
        for (size_t i = 0; i < res.atoms.size(); ++i) {
            EXPECT_NEAR(coords[i].x, before_coords[i][0], 1e-3);
            EXPECT_NEAR(coords[i].y, before_coords[i][1], 1e-3);
            EXPECT_NEAR(coords[i].z, before_coords[i][2], 1e-3);
        }
        
        // Apply rotation of 30 degrees
        // Rotation about bond: CA (index 1) -> CB (index 4)
        // Rotating subtree of CB (index 4)
        int result = coot::torsion_general::change_by_gemmi(
            30.0, coords, atom_CA_idx, atom_CB_idx, children_of_CB);
        
        EXPECT_EQ(result, 0);
        
        // AFTER coordinates (from oracle - after 30 degree rotation)
        double after_coords[12][3] = {
            {17.343, 12.361, 58.028},        // N (unchanged)
            {16.400, 11.464, 58.677},        // CA (unchanged)
            {16.184, 11.852, 60.139},        // C (unchanged)
            {15.067, 12.195, 60.529},        // O (unchanged)
            {16.830, 10.003, 58.517},        // CB (unchanged)
            {17.1579, 9.59855, 57.1006},     // CG
            {16.2407, 9.78164, 56.0737},     // CD1
            {18.3801, 9.02308, 56.7856},     // CD2
            {16.531, 9.41198, 54.7695},      // CE1
            {18.6861, 8.64526, 55.4866},     // CE2
            {17.7578, 8.84016, 54.476},      // CZ
            {18.0555, 8.46898, 53.1873},     // OH
        };
        
        for (size_t i = 0; i < res.atoms.size(); ++i) {
            EXPECT_NEAR(coords[i].x, after_coords[i][0], 1e-3);
            EXPECT_NEAR(coords[i].y, after_coords[i][1], 1e-3);
            EXPECT_NEAR(coords[i].z, after_coords[i][2], 1e-3);
        }
    }
    
    // Case 2: diff=0 — should still return 0 but coordinates unchanged
    {
        // Reset coordinates to original (unrotated)
        coords.clear();
        for (const auto& atom : res.atoms) {
            coords.push_back(atom.pos);
        }
        
        // BEFORE coordinates (same as original PDB)
        double before_coords[12][3] = {
            {17.343, 12.361, 58.028},
            {16.400, 11.464, 58.677},
            {16.184, 11.852, 60.139},
            {15.067, 12.195, 60.529},
            {16.830, 10.003, 58.517},
            {16.483, 9.391, 57.182},
            {15.177, 9.399, 56.709},
            {17.456, 8.794, 56.394},
            {14.847, 8.838, 55.485},
            {17.143, 8.226, 55.168},
            {15.834, 8.248, 54.713},
            {15.518, 7.688, 53.499},
        };
        
        // Verify BEFORE coordinates match PDB
        for (size_t i = 0; i < res.atoms.size(); ++i) {
            EXPECT_NEAR(coords[i].x, before_coords[i][0], 1e-3);
            EXPECT_NEAR(coords[i].y, before_coords[i][1], 1e-3);
            EXPECT_NEAR(coords[i].z, before_coords[i][2], 1e-3);
        }
        
        // Apply rotation of 0 degrees
        int result = coot::torsion_general::change_by_gemmi(
            0.0, coords, atom_CA_idx, atom_CB_idx, children_of_CB);
        
        EXPECT_EQ(result, 0);
        
        // AFTER coordinates should be identical to BEFORE (unchanged)
        for (size_t i = 0; i < res.atoms.size(); ++i) {
            EXPECT_NEAR(coords[i].x, before_coords[i][0], 1e-3);
            EXPECT_NEAR(coords[i].y, before_coords[i][1], 1e-3);
            EXPECT_NEAR(coords[i].z, before_coords[i][2], 1e-3);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}