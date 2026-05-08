#include "function.hh"
#include <gemmi/neighbor.hpp>
#include <gemmi/calculate.hpp>

namespace coot {

void molecule_t::add_neighbor_residues_for_refinement_help_gemmi(gemmi::Structure& st) {
    neighbouring_residues.clear();
    
    if (st.models.empty() || st.models[0].chains.empty()) {
        return;
    }
    
    // Build neighbor search
    gemmi::NeighborSearch ns(st.models[0], st.cell, 5.0);
    ns.populate(false);
    
    // Collect all residues from the structure
    std::vector<std::pair<bool, gemmi::Residue*>> all_residues;
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                all_residues.push_back({true, &residue});
            }
        }
    }
    
    if (all_residues.empty()) {
        return;
    }
    
    // Find neighboring residues (same logic as original: residues within 5Å of any residue in the fragment)
    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> rnr;
    
    for (size_t i = 0; i < all_residues.size(); ++i) {
        gemmi::Residue* res = all_residues[i].second;
        std::set<gemmi::Residue*> neighbors;
        
        // Find atoms in this residue
        for (auto& atom : res->atoms) {
            auto hits = ns.find_atoms(atom.pos, '\0', 0.0, 5.0);
            for (auto* mark : hits) {
                auto cra = mark->to_cra(st.models[0]);
                if (cra.residue && cra.residue != res) {
                    neighbors.insert(cra.residue);
                }
            }
        }
        
        if (!neighbors.empty()) {
            rnr[res] = neighbors;
        }
    }
    
    // Convert map of sets to vector of pairs (keeping original logic)
    std::set<gemmi::Residue*> keys;
    for (const auto& kv : rnr) {
        keys.insert(kv.first);
    }
    
    for (const auto& kv : rnr) {
        for (gemmi::Residue* neighbor : kv.second) {
            if (keys.find(neighbor) == keys.end()) {
                // Not in the original set, so it's a neighboring residue
                neighbouring_residues.push_back({false, neighbor});
            }
        }
    }
}

}