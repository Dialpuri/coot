#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <vector>
#include <cmath>
#include <iostream>

TEST(OracleTest, AtomTreeConstructor) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    std::vector<gemmi::Atom> atoms;
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name == "A") {
                for (const auto& res : chain.residues) {
                    for (const auto& atom : res.atoms) {
                        atoms.push_back(atom);
                    }
                }
            }
        }
    }

    std::cout << "INPUT n_selected_atoms: " << atoms.size() << std::endl;

    std::vector<std::vector<int>> contact_indices;
    contact_indices.resize(atoms.size());
    
    for (size_t i = 0; i < atoms.size() - 1; i++) {
        if (atoms[i].name[0] == 'C' && atoms[i].name[1] == 'A') {
            for (size_t j = i + 1; j < atoms.size(); j++) {
                if (atoms[j].name[0] == 'C' && atoms[j].name[1] == 'A') {
                    double dx = atoms[i].pos.x - atoms[j].pos.x;
                    double dy = atoms[i].pos.y - atoms[j].pos.y;
                    double dz = atoms[i].pos.z - atoms[j].pos.z;
                    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (dist < 4.0) {
                        contact_indices[i].push_back(j);
                        contact_indices[j].push_back(i);
                        break;
                    }
                }
            }
        }
    }

    {
        int base_atom_index = 0;
        for (size_t i = 0; i < atoms.size(); i++) {
            if (atoms[i].name[0] == 'C' && atoms[i].name[1] == 'A') {
                base_atom_index = i;
                break;
            }
        }

        std::cout << "INPUT base_atom_index: " << base_atom_index << std::endl;
        std::cout << "INPUT contact_indices_size: " << contact_indices.size() << std::endl;

        coot::atom_tree_t tree(contact_indices, base_atom_index, atoms, 0);
        
        std::cout << "OUTPUT tree_constructed: true" << std::endl;
        std::cout << "OUTPUT atom_vertex_vec_size: " << tree.atom_vertex_vec.size() << std::endl;
    }

    {
        std::vector<std::vector<int>> empty_contacts;
        int base_atom_index = 0;
        
        std::cout << "INPUT empty_contacts: true" << std::endl;
        std::cout << "INPUT base_atom_index: " << base_atom_index << std::endl;

        coot::atom_tree_t tree(empty_contacts, base_atom_index, atoms, 0);
        
        std::cout << "OUTPUT empty_test_constructed: true" << std::endl;
        std::cout << "OUTPUT empty_test_atom_vertex_vec_size: " << tree.atom_vertex_vec.size() << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}