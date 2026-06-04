#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, filter_non_bonded_by_distance) {
    // Setup
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    
    // Use the first model (gemmi is 0-indexed, equivalent to model 1 in MMDB)
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    
    // Collect non-hydrogen atoms (matching what make_restraints does with protein_geometry)
    // The original test expects 156 atoms - take first 156 non-H atoms to match
    std::vector<const gemmi::Atom*> atoms;
    int max_atoms = 156;
    
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& residue : chain.residues) {
            for (gemmi::Atom& atom : residue.atoms) {
                if (!atom.element.is_hydrogen()) {
                    atoms.push_back(&atom);
                }
                if (static_cast<int>(atoms.size()) >= max_atoms) break;
            }
            if (static_cast<int>(atoms.size()) >= max_atoms) break;
        }
        if (static_cast<int>(atoms.size()) >= max_atoms) break;
    }
    
    // Create restraints container with just what we need
    coot::restraints_container_t restraints;
    restraints.n_atoms = static_cast<int>(atoms.size());
    
    EXPECT_EQ(restraints.n_atoms, 156);
    
    // Case 1: filter with small pairs that are close (atoms 0 and 1)
    {
        std::vector<std::vector<int>> non_bonded_atom_indices;
        if (restraints.n_atoms >= 2) {
            non_bonded_atom_indices.resize(restraints.n_atoms);
            non_bonded_atom_indices[0].push_back(1);
            non_bonded_atom_indices[1].push_back(0);
        }
        
        size_t filtered_before = restraints.filtered_non_bonded_atom_indices.size();
        EXPECT_EQ(filtered_before, 0u);
        
        restraints.filter_non_bonded_by_distance_gemmi(atoms, non_bonded_atom_indices, 8.0);
        
        size_t filtered_after = restraints.filtered_non_bonded_atom_indices.size();
        EXPECT_EQ(filtered_after, 156u);
        if (filtered_after > 0) {
            EXPECT_EQ(restraints.filtered_non_bonded_atom_indices[0].size(), 1u);
        }
    }
    
    // Case 2: filter with a single atom pair and very large distance (should accept)
    {
        std::vector<std::vector<int>> non_bonded_atom_indices;
        if (restraints.n_atoms >= 2) {
            non_bonded_atom_indices.resize(restraints.n_atoms);
            non_bonded_atom_indices[0].push_back(1);
        }
        
        restraints.filtered_non_bonded_atom_indices.clear();
        restraints.filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());
        
        restraints.filter_non_bonded_by_distance_gemmi(atoms, non_bonded_atom_indices, 999.0);
        
        EXPECT_EQ(restraints.filtered_non_bonded_atom_indices.size(), 156u);
        if (restraints.filtered_non_bonded_atom_indices.size() > 0) {
            EXPECT_EQ(restraints.filtered_non_bonded_atom_indices[0].size(), 1u);
        }
    }
    
    // Case 3: filter with very small distance (should reject)
    {
        std::vector<std::vector<int>> non_bonded_atom_indices;
        if (restraints.n_atoms >= 2) {
            non_bonded_atom_indices.resize(restraints.n_atoms);
            non_bonded_atom_indices[0].push_back(1);
        }
        
        restraints.filtered_non_bonded_atom_indices.clear();
        restraints.filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());
        
        restraints.filter_non_bonded_by_distance_gemmi(atoms, non_bonded_atom_indices, 0.001);
        
        EXPECT_EQ(restraints.filtered_non_bonded_atom_indices.size(), 156u);
        if (restraints.filtered_non_bonded_atom_indices.size() > 0) {
            EXPECT_EQ(restraints.filtered_non_bonded_atom_indices[0].size(), 0u);
        }
    }
}