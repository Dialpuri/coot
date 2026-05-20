#include "function.hh"
#include <algorithm>
#include <iostream>

namespace coot {

void match_container_for_residues_t::delete_downstream_gemmi(
    gemmi::Structure& st, bool from_first,
    const std::vector<gemmi::CRA>& atoms0,
    const std::vector<gemmi::CRA>& atoms1) {

    bool debug = false;

    // Select the appropriate atom vector based on from_first
    const std::vector<gemmi::CRA>& atoms = from_first ? atoms0 : atoms1;
    
    if (atoms.empty()) return;
    
    // Get the first atom's chain to identify which chain we're working on
    const gemmi::Chain* chain_p = atoms[0].chain;
    if (!chain_p) return;
    
    // Find matching residue
    gemmi::Residue* matchers_residue = nullptr;
    bool found_matchers = false;
    
    // Check if any atom is part of a matching pair
    for (const auto& cra : atoms) {
        if (from_first) {
            for (size_t ip = 0; ip < atom_pairs.size(); ip++) {
                if (atom_pairs[ip].first.chain == cra.chain &&
                    atom_pairs[ip].first.residue == cra.residue &&
                    atom_pairs[ip].first.atom == cra.atom) {
                    found_matchers = true;
                    matchers_residue = cra.residue;
                    break;
                }
            }
        } else {
            for (size_t ip = 0; ip < atom_pairs.size(); ip++) {
                if (atom_pairs[ip].second.chain == cra.chain &&
                    atom_pairs[ip].second.residue == cra.residue &&
                    atom_pairs[ip].second.atom == cra.atom) {
                    found_matchers = true;
                    matchers_residue = cra.residue;
                    break;
                }
            }
        }
        if (found_matchers) break;
    }
    
    if (!found_matchers) return;
    
    // Collect residues to delete (downstream from matchers_residue)
    std::vector<gemmi::Residue*> delete_these_residues;
    
    for (const auto& cra : atoms) {
        if (cra.chain != chain_p) continue;  // Only same chain
        
        if (cra.residue != matchers_residue) {
            // Check if already in delete list
            bool already = false;
            for (auto* r : delete_these_residues) {
                if (r == cra.residue) {
                    already = true;
                    break;
                }
            }
            if (!already) {
                delete_these_residues.push_back(cra.residue);
            }
        }
    }
    
    // Delete the residues
    if (delete_these_residues.size() > 0) {
        if (debug) {
            for (auto* r : delete_these_residues) {
                std::cout << "debug:: in delete_downstream_gemmi() delete "
                          << r->name << " " << r->seqid.num.value << std::endl;
            }
        }
        
        // Remove residues from chain (cast away const since we need to modify)
        std::vector<gemmi::Residue>& residues = const_cast<std::vector<gemmi::Residue>&>(chain_p->residues);
        for (auto* r : delete_these_residues) {
            auto it = std::find(residues.begin(), residues.end(), *r);
            if (it != residues.end()) {
                if (debug) {
                    std::cout << "debug:: in delete_downstream_gemmi() about to delete residue "
                              << r->name << " " << r->seqid.num.value << std::endl;
                }
                residues.erase(it);
            }
        }
    }
}

} // namespace coot