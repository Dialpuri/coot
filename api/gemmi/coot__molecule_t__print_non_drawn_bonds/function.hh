#pragma once
#include <set>
#include <string>
#include <iostream>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/calculate.hpp>

namespace coot {

class molecule_t {
public:
    std::set<int> no_bonds_to_these_atom_indices;
    gemmi::Structure structure;

    // Constructor that loads a structure from a PDB file
    molecule_t(const std::string& path) : structure(gemmi::read_pdb_file(path)) {}

    // Print atoms that have no bonds drawn
    void print_non_drawn_bonds() const {
        std::cout << "----------- no bonds to these atoms table: " << std::endl;
        size_t total = gemmi::count_atom_sites(structure);
        for (int idx : no_bonds_to_these_atom_indices) {
            if (idx >= 0) {
                if (static_cast<size_t>(idx) < total) {
                    // Find the atom at this index
                    size_t current = 0;
                    for (const gemmi::Model& model : structure.models) {
                        for (const gemmi::Chain& chain : model.chains) {
                            for (const gemmi::Residue& res : chain.residues) {
                                for (const gemmi::Atom& atom : res.atoms) {
                                    if (current == static_cast<size_t>(idx)) {
                                        std::cout << "  " << idx << "  "
                                                  << chain.name << " "
                                                  << res.seqid.num.value << res.seqid.icode << " "
                                                  << atom.name << " "
                                                  << std::string(1, atom.altloc) << std::endl;
                                        goto next_idx;
                                    }
                                    current++;
                                }
                            }
                        }
                    }
                next_idx: {}
                } else {
                    std::cout << "ERROR:: atom index " << idx << " out of range" << std::endl;
                }
            }
        }
    }
};

} // namespace coot