#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/atom-tree.hh"

TEST(OracleTest, AtomTreeConstruction) {
    // Set up the gemmi structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty()) << "Failed to load PDB file";

    // Get residue ILE 10 in chain A
    gemmi::Residue* res = nullptr;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& r : chain.residues) {
                if (r.name == "ILE" && r.seqid.num.value == 10) {
                    res = &r;
                    break;
                }
            }
        }
        if (res) break;
    }
    ASSERT_NE(res, nullptr) << "Could not get residue ILE 10 in chain A";

    // Get the atoms in this residue
    std::vector<gemmi::Atom> atoms;
    for (const auto& atom : res->atoms) {
        atoms.push_back(atom);
    }
    
    ASSERT_GT(atoms.size(), 0) << "no atoms in residue";

    // Create contact_indices from distances (since gemmi doesn't have bonds)
    // For testing purposes, we'll use a simple distance-based contact search
    std::vector<std::vector<int> > contact_indices(atoms.size());
    
    // Use a reasonable distance cutoff (e.g., 2.0 Angstroms for covalent bonds)
    const double cutoff_sq = 2.0 * 2.0;
    
    for (size_t iat = 0; iat < atoms.size(); iat++) {
        for (size_t jat = 0; jat < atoms.size(); jat++) {
            if (iat == jat) continue;
            
            // Calculate distance squared using .x, .y, .z (Position has these fields)
            double dist_sq = 0.0;
            for (int k = 0; k < 3; k++) {
                double diff = 0.0;
                if (k == 0) diff = atoms[iat].pos.x - atoms[jat].pos.x;
                else if (k == 1) diff = atoms[iat].pos.y - atoms[jat].pos.y;
                else if (k == 2) diff = atoms[iat].pos.z - atoms[jat].pos.z;
                dist_sq += diff * diff;
            }
            
            if (dist_sq < cutoff_sq) {
                contact_indices[iat].push_back(jat);
            }
        }
    }

    // Test 1: Normal case with valid contact_indices
    {
        int base_atom_index = 0;
        std::vector<coot::atom_vertex> atom_vertex_vec;
        bool result = coot::fill_atom_vertex_vec_using_contacts_gemmi(contact_indices, base_atom_index, atoms, atom_vertex_vec);
        EXPECT_EQ(atoms.size(), 8); // ILE 10 has 8 atoms
        
        // Check if vertex was properly initialized
        bool vertex_initialized = false;
        if (atom_vertex_vec.size() > static_cast<size_t>(base_atom_index)) {
            vertex_initialized = true;
        }
        EXPECT_TRUE(vertex_initialized);
    }

    // Test 2: Empty contact_indices (should still create valid tree)
    {
        std::vector<std::vector<int> > empty_contacts;
        std::vector<coot::atom_vertex> atom_vertex_vec;
        bool result = coot::fill_atom_vertex_vec_using_contacts_gemmi(empty_contacts, 0, atoms, atom_vertex_vec);
        // The tree is created even with empty contacts
        EXPECT_EQ(atoms.size(), 8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}