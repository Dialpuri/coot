#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Test fixture using the example.pdb file
TEST(OracleTest, NewPositionsForResidueAtomsGemmi) {
    // Load the test PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";
    
    gemmi::Model& model = st.models[0];
    
    // Test 1: Valid residue move with CID //A/10
    std::string residue_cid = "//A/10";
    
    // In gemmi, we need to find the residue by chain name and sequence number
    // CID format is //chain/resnum - parse it
    // For "//A/10", we need chain "A" and residue number 10
    gemmi::Chain* chain = nullptr;
    for (auto& c : model.chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr) << "Failed to find chain A";
    
    // Find residue with seqid.num.value == 10
    gemmi::Residue* res = nullptr;
    for (auto& r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res = &r;
            break;
        }
    }
    ASSERT_NE(res, nullptr) << "Failed to find residue with seqid 10 in chain A";
    
    // Get original positions BEFORE move
    int natoms = res->atoms.size();
    ASSERT_GT(natoms, 0) << "Residue has no atoms";
    
    std::vector<float> original_x, original_y, original_z;
    for (int i = 0; i < std::min(3, natoms); i++) {
        original_x.push_back(res->atoms[i].pos.x);
        original_y.push_back(res->atoms[i].pos.y);
        original_z.push_back(res->atoms[i].pos.z);
    }
    
    // Create moved atoms with new positions (offset by +1 in x direction)
    // In gemmi, we need to map atom names to atoms
    std::vector<coot::api::moved_atom_t> moved_atoms;
    
    for (int i = 0; i < natoms; i++) {
        const auto& atom = res->atoms[i];
        std::string atom_name = atom.name;
        
        // New position - move by 1 Angstrom in x direction
        float new_x = atom.pos.x + 1.0f;
        float new_y = atom.pos.y + 0.0f;
        float new_z = atom.pos.z + 0.0f;
        
        coot::api::moved_atom_t ma(atom_name, "", new_x, new_y, new_z);
        moved_atoms.push_back(ma);
    }
    
    int status = coot::molecule_t::new_positions_for_residue_atoms_gemmi(
        model, residue_cid, moved_atoms);
    
    // Status should be positive (number of atoms processed)
    EXPECT_EQ(status, natoms) << "Status should equal number of atoms moved";
    
    // Verify atoms were moved by checking positions again
    // Re-find the residue (in case it was modified)
    gemmi::Residue* res2 = nullptr;
    for (auto& r : chain->residues) {
        if (r.seqid.num.value == 10) {
            res2 = &r;
            break;
        }
    }
    ASSERT_NE(res2, nullptr) << "Residue disappeared after move";
    
    // Check first 3 atoms have been moved by +1 in x direction
    for (int i = 0; i < std::min(3, natoms); i++) {
        const auto& atom = res2->atoms[i];
        EXPECT_FLOAT_EQ(atom.pos.x, original_x[i] + 1.0f) 
            << "Atom " << atom.name << " x-coordinate should be moved by 1.0";
        EXPECT_FLOAT_EQ(atom.pos.y, original_y[i]) 
            << "Atom " << atom.name << " y-coordinate should be unchanged";
        EXPECT_FLOAT_EQ(atom.pos.z, original_z[i]) 
            << "Atom " << atom.name << " z-coordinate should be unchanged";
    }
    
    // Test 2: Invalid CID to test error case
    std::string invalid_cid = "//Z/999";
    
    std::vector<coot::api::moved_atom_t> empty_moved;
    int status2 = coot::molecule_t::new_positions_for_residue_atoms_gemmi(
        model, invalid_cid, empty_moved);
    
    EXPECT_EQ(status2, -1) << "Invalid CID should return -1";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}