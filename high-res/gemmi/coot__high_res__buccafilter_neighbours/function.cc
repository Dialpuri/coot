#include "function.hh"
#include <gemmi/polyheur.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <algorithm>

namespace coot {
namespace high_res {

static void mark_neighbours_gemmi(int iatom, int igroup, const std::string &atom_name,
                           const std::vector<std::vector<int> > &neighbours,
                           const std::vector<gemmi::Atom*> &atom_selection,
                           std::map<gemmi::Atom*, int> &group_map) {
    // Mark the current atom
    if (iatom < 0 || iatom >= static_cast<int>(atom_selection.size())) {
        return;
    }
    gemmi::Atom* atom = atom_selection[iatom];
    if (group_map[atom] != -1) {
        return; // already marked
    }
    group_map[atom] = igroup;
    
    // Recursively mark all neighbors with the same atom name
    if (iatom < static_cast<int>(neighbours.size())) {
        for (int neighbor_idx : neighbours[iatom]) {
            if (neighbor_idx >= 0 && neighbor_idx < static_cast<int>(atom_selection.size())) {
                gemmi::Atom* neighbor = atom_selection[neighbor_idx];
                if (neighbor->name == atom_name && group_map[neighbor] == -1) {
                    mark_neighbours_gemmi(neighbor_idx, igroup, atom_name, 
                                         neighbours, atom_selection, group_map);
                }
            }
        }
    }
}

void buccafilter_neighbours_gemmi(gemmi::Structure& st) {
    // Ensure entities are set up
    gemmi::setup_entities(st);
    
    if (st.models.empty()) {
        return;
    }
    
    gemmi::Model& model = st.models[0];
    
    // Collect all atoms into a vector for processing
    std::vector<gemmi::Atom*> atom_selection;
    for (auto& chain : model.chains) {
        for (auto& res : chain.residues) {
            for (auto& atom : res.atoms) {
                atom_selection.push_back(&atom);
            }
        }
    }
    
    if (atom_selection.empty()) {
        return;
    }
    
    // Find contacts within 2.0 Å using NeighborSearch
    const float max_dist = 2.0f;
    gemmi::NeighborSearch ns(model, st.cell, max_dist);
    ns.populate(false); // false = don't include hydrogens
    
    // Build neighbor list as vector of vectors of indices
    std::vector<std::vector<int>> neighbours(atom_selection.size());
    
    for (size_t i = 0; i < atom_selection.size(); ++i) {
        const gemmi::Position pos = atom_selection[i]->pos;
        std::vector<gemmi::NeighborSearch::Mark*> marks = 
            ns.find_atoms(pos, '\0', 0.0, max_dist);
        
        for (auto* mark : marks) {
            gemmi::CRA cra = mark->to_cra(model);
            if (cra.atom && cra.atom != atom_selection[i]) {
                // Find the index of the neighbor
                for (size_t j = 0; j < atom_selection.size(); ++j) {
                    if (atom_selection[j] == cra.atom) {
                        neighbours[i].push_back(static_cast<int>(j));
                        break;
                    }
                }
            }
        }
    }
    
    // Assign group IDs to atoms using UDD-like integer storage
    std::map<gemmi::Atom*, int> group_map;
    for (auto* atom : atom_selection) {
        group_map[atom] = -1; // -1 means unset
    }
    
    // Group atoms by connectivity (same atom name → neighbors → same group)
    int igroup = 0;
    for (size_t iat = 0; iat < atom_selection.size(); ++iat) {
        if (group_map[atom_selection[iat]] == -1) {
            // Mark this atom and all connected atoms with same name
            std::string atom_name = atom_selection[iat]->name;
            mark_neighbours_gemmi(static_cast<int>(iat), igroup, 
                                  atom_name, neighbours, atom_selection, group_map);
            igroup++;
        }
    }
    
    // Build groups from the group_map
    std::vector<std::vector<int>> groups(igroup);
    for (size_t iat = 0; iat < atom_selection.size(); ++iat) {
        int gid = group_map[atom_selection[iat]];
        if (gid >= 0 && gid < igroup) {
            groups[gid].push_back(static_cast<int>(iat));
        }
    }
    
    // Filter on groups - keep only the largest group
    if (!groups.empty()) {
        // Find the largest group
        size_t largest_idx = 0;
        size_t largest_size = 0;
        for (size_t i = 0; i < groups.size(); ++i) {
            if (groups[i].size() > largest_size) {
                largest_size = groups[i].size();
                largest_idx = i;
            }
        }
        
        // Keep only atoms in the largest group
        std::vector<bool> keep(atom_selection.size(), false);
        for (int idx : groups[largest_idx]) {
            keep[idx] = true;
        }
        
        // Remove atoms not in the largest group
        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                size_t write_pos = 0;
                for (size_t read_pos = 0; read_pos < res.atoms.size(); ++read_pos) {
                    gemmi::Atom* atom = &res.atoms[read_pos];
                    // Find index in atom_selection
                    int idx = -1;
                    for (size_t i = 0; i < atom_selection.size(); ++i) {
                        if (atom_selection[i] == atom) {
                            idx = static_cast<int>(i);
                            break;
                        }
                    }
                    if (idx >= 0 && keep[static_cast<size_t>(idx)]) {
                        if (read_pos != write_pos) {
                            res.atoms[write_pos] = std::move(res.atoms[read_pos]);
                        }
                        write_pos++;
                    }
                }
                res.atoms.resize(write_pos);
            }
        }
        
        // Remove empty residues and chains
        for (auto& chain : model.chains) {
            chain.residues.erase(
                std::remove_if(chain.residues.begin(), chain.residues.end(),
                    [](const gemmi::Residue& res) { return res.atoms.empty(); }),
                chain.residues.end()
            );
        }
        model.chains.erase(
            std::remove_if(model.chains.begin(), model.chains.end(),
                [](const gemmi::Chain& chain) { return chain.residues.empty(); }),
            model.chains.end()
        );
    }
}

} // namespace high_res
} // namespace coot