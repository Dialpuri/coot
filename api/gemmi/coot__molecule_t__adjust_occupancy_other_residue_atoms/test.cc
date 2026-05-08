#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, adjust_occupancy_other_residue_atoms) {
    // Load PDB using gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get first model, first chain, first residue
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& residue = chain.residues[0];
    
    // Get original atom (e.g., N)
    gemmi::Atom& original_atom = residue.atoms[0];
    char original_altloc = original_atom.altloc;
    
    // Create a new atom with same name but different altLoc
    gemmi::Atom new_atom = original_atom;
    new_atom.name = original_atom.name;
    new_atom.altloc = 'A';  // set altLoc to 'A'
    new_atom.pos.x += 0.5;
    new_atom.pos.y += 0.5;
    new_atom.pos.z += 0.5;
    new_atom.occ = 0.6f;  // >0.5 to trigger adjustment
    residue.atoms.push_back(new_atom);
    
    // Get pointer to the newly added atom
    gemmi::Atom* new_atom_ptr = &residue.atoms.back();
    
    // Call the function with the new atom
    coot::molecule_t::adjust_occupancy_other_residue_atoms_gemmi(
        new_atom_ptr, &residue, false);
    
    // Record all atoms after
    fprintf(stderr, "All atoms after function call:\n");
    std::vector<std::pair<float,char>> atom_occ_after;
    for (const auto& atom : residue.atoms) {
        fprintf(stderr, "  name='%s', altloc='%c', occ=%.3f\n", 
                atom.name.c_str(), atom.altloc, atom.occ);
        if (atom.name == original_atom.name) {
            atom_occ_after.push_back({atom.occ, atom.altloc});
        }
    }
    
    // Check that we have 2 atoms with same name (original + new)
    EXPECT_EQ(atom_occ_after.size(), 2) << "Should have 2 atoms with same name";
    
    // Verify that at least one atom's occupancy was modified
    bool modified = false;
    if (atom_occ_after.size() == 2) {
        // Check if any atom has occupancy 0.4
        for (const auto& p : atom_occ_after) {
            if (std::abs(p.first - 0.4f) < 0.001f) {
                modified = true;
                break;
            }
        }
    }
    EXPECT_TRUE(modified) << "Atom occupancies were not modified correctly";
    
    // Check that there is an atom with blank altLoc that was adjusted to 0.4
    bool found_blank_altloc_with_0_4 = false;
    bool found_target_atom_with_0_6 = false;
    for (const auto& p : atom_occ_after) {
        if ((p.second == ' ' || p.second == '\0') && std::abs(p.first - 0.4f) < 0.001f) {
            found_blank_altloc_with_0_4 = true;
        }
        if (p.second == 'A' && std::abs(p.first - 0.6f) < 0.001f) {
            found_target_atom_with_0_6 = true;
        }
    }
    
    EXPECT_TRUE(found_blank_altloc_with_0_4) << "Expected an atom with blank altLoc and occupancy ~0.4";
    EXPECT_TRUE(found_target_atom_with_0_6) << "Expected a target atom with altLoc='A' and occupancy ~0.6";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}