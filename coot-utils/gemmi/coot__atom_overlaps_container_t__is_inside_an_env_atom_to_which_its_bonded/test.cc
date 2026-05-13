#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <clipper/core/coords.h>
#include "function.hh"

// Test for atom_overlaps_container_t::is_inside_an_env_atom_to_which_its_bonded_gemmi
TEST(AtomOverlapsTest, IsInsideAnEnvAtomToWhichItsBonded) {
    // Load a simple PDB file for testing
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Get the first model
    gemmi::Model& model = st.models[0];
    
    // Collect atoms from a specific residue for testing
    std::vector<gemmi::Atom*> env_residue_atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "ILE" && res.seqid.num.value == 10) {
                for (gemmi::Atom& atom : res.atoms) {
                    env_residue_atoms.push_back(&atom);
                }
                break;
            }
        }
    }
    
    // Create an atom_overlaps_container_t instance
    coot::atom_overlaps_container_t aoc;
    aoc.neighb_atom_radius = {1.6, 1.6, 1.6}; // dummy radii for testing
    
    // Test with a point inside an atom
    // Get coordinates of first atom
    clipper::Coord_orth pt_at_surface(env_residue_atoms[0]->pos.x,
                                      env_residue_atoms[0]->pos.y,
                                      env_residue_atoms[0]->pos.z);
    
    // Test bonded neighbors (indices 0, 1, 2)
    std::vector<int> bonded_neighb_indices = {0, 1, 2};
    
    // The point should be "inside" the first atom (distance = 0)
    bool result = aoc.is_inside_an_env_atom_to_which_its_bonded_gemmi(
        0, bonded_neighb_indices, env_residue_atoms, pt_at_surface);
    
    // Since the point is exactly at the atom center, it should be considered "inside"
    EXPECT_TRUE(result);
}

// Test with a point outside all atoms
TEST(AtomOverlapsTest, IsOutsideAllEnvAtoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    
    std::vector<gemmi::Atom*> env_residue_atoms;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "ILE" && res.seqid.num.value == 10) {
                for (gemmi::Atom& atom : res.atoms) {
                    env_residue_atoms.push_back(&atom);
                }
                break;
            }
        }
    }
    
    coot::atom_overlaps_container_t aoc;
    aoc.neighb_atom_radius = {1.6, 1.6, 1.6};
    
    // Create a point far away from all atoms (more than 5 Angstroms from any atom)
    clipper::Coord_orth pt_far_away(100.0, 100.0, 100.0);
    std::vector<int> bonded_neighb_indices = {0, 1, 2};
    
    bool result = aoc.is_inside_an_env_atom_to_which_its_bonded_gemmi(
        0, bonded_neighb_indices, env_residue_atoms, pt_far_away);
    
    EXPECT_FALSE(result);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}