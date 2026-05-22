#include <gemmi/neighbor.hpp>
#include <vector>
#include <string>
#include "function.hh"

namespace coot {

std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>>
interface_residues_gemmi(gemmi::Structure& st,
                          const std::string& chain_A,
                          const std::string& chain_B,
                          float min_dist) {
    std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>> result;
    
    if (st.models.empty())
        return result;
    
    gemmi::Model& model = st.models[0];
    
    // Build neighbor search
    gemmi::NeighborSearch ns(model, st.cell, min_dist);
    ns.populate(false);
    
    // Get residues for each chain
    std::vector<gemmi::Residue*> residues_A;
    std::vector<gemmi::Residue*> residues_B;
    
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_A) {
            for (gemmi::Residue& res : chain.residues) {
                residues_A.push_back(&res);
            }
        }
        if (chain.name == chain_B) {
            for (gemmi::Residue& res : chain.residues) {
                residues_B.push_back(&res);
            }
        }
    }
    
    // Build all_molecule_map: only inter-chain contacts
    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> all_molecule_map;
    
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            for (const auto& atom : res.atoms) {
                auto hits = ns.find_atoms(atom.pos, '\0', 0.0, min_dist);
                for (auto* mark : hits) {
                    gemmi::CRA cra = mark->to_cra(model);
                    if (cra.residue && cra.residue != &res) {
                        // Only inter-chain contacts (different chain names)
                        if (cra.chain->name != chain.name) {
                            all_molecule_map[&res].insert(cra.residue);
                        }
                    }
                }
            }
        }
    }
    
    // Filter: keep only entries where key is in the limit set
    auto filter_by_limits = [](const std::map<gemmi::Residue*, std::set<gemmi::Residue*>>& m,
                               const std::set<gemmi::Residue*>& limit) {
        std::map<gemmi::Residue*, std::set<gemmi::Residue*>> result;
        for (const auto& [key, neighbors] : m) {
            if (limit.count(key)) {
                result[key] = neighbors;
            }
        }
        return result;
    };
    
    std::set<gemmi::Residue*> limit_A(residues_A.begin(), residues_A.end());
    std::set<gemmi::Residue*> limit_B(residues_B.begin(), residues_B.end());
    
    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> mA = filter_by_limits(all_molecule_map, limit_A);
    std::map<gemmi::Residue*, std::set<gemmi::Residue*>> mB = filter_by_limits(all_molecule_map, limit_B);
    
    // Extract keys
    std::set<gemmi::Residue*> r_1_residues;
    std::set<gemmi::Residue*> r_2_residues;
    
    for (const auto& [key, _] : mA) {
        r_1_residues.insert(key);
    }
    for (const auto& [key, _] : mB) {
        r_2_residues.insert(key);
    }
    
    result = std::make_pair(r_1_residues, r_2_residues);
    return result;
}

} // namespace coot