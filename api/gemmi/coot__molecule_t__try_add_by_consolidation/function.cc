#include "function.hh"
#include <gemmi/pdb.hpp>
#include <map>
#include <algorithm>
#include <set>

namespace coot {

std::pair<bool, std::vector<std::string>>
molecule_t_try_add_by_consolidation_gemmi(const gemmi::Structure& source_structure,
                                          gemmi::Structure& target_structure) {
    bool status = false;
    std::vector<std::string> chain_ids;

    // Build a map of residue_type -> chain_id for target_structure
    // For each model and chain in target, collect chains with single residue type
    std::map<std::string, std::string> single_res_type_map;
    
    for (const gemmi::Model& model : target_structure.models) {
        for (const gemmi::Chain& chain : model.chains) {
            std::vector<std::string> residue_types;
            for (const gemmi::Residue& res : chain.residues) {
                std::string res_name = res.name;
                if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end()) {
                    residue_types.push_back(res_name);
                }
                if (residue_types.size() > 1) break;
            }
            if (residue_types.size() == 1) {
                single_res_type_map[residue_types[0]] = chain.name;
            }
        }
    }

    // Process each model and chain in source_structure
    for (const gemmi::Model& source_model : source_structure.models) {
        for (const gemmi::Chain& source_chain : source_model.chains) {
            bool done_this_chain = false;
            std::vector<std::string> residue_types;
            for (const gemmi::Residue& res : source_chain.residues) {
                std::string res_name = res.name;
                if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end()) {
                    residue_types.push_back(res_name);
                }
            }

            if (residue_types.size() == 1 && residue_types[0] != "ALA") {
                auto it = single_res_type_map.find(residue_types[0]);
                if (it != single_res_type_map.end()) {
                    // Match found - copy chain from source to target
                    std::string target_chain_id = it->second;
                    
                    // Add residues from source chain to target chain
                    for (gemmi::Model& model : target_structure.models) {
                        for (gemmi::Chain& chain : model.chains) {
                            if (chain.name == target_chain_id) {
                                for (const gemmi::Residue& res : source_chain.residues) {
                                    chain.residues.push_back(res);
                                }
                                done_this_chain = true;
                                if (std::find(chain_ids.begin(), chain_ids.end(), target_chain_id) == chain_ids.end()) {
                                    chain_ids.push_back(target_chain_id);
                                }
                            }
                        }
                        if (done_this_chain) break;
                    }
                }
            }

            if (!done_this_chain) {
                // Copy whole chain to a new chain with new ID
                std::string current_chain_id = source_chain.name;
                std::string new_chain_id = current_chain_id;
                
                // Generate a new chain ID if it already exists in target
                std::set<std::string> existing_chain_ids;
                for (const gemmi::Model& model : target_structure.models) {
                    for (const gemmi::Chain& chain : model.chains) {
                        existing_chain_ids.insert(chain.name);
                    }
                }
                
                // Try to find a unique chain ID
                if (existing_chain_ids.count(new_chain_id) > 0) {
                    // Simple approach: append a number
                    for (int i = 1; i < 100; ++i) {
                        new_chain_id = current_chain_id + std::to_string(i);
                        if (existing_chain_ids.count(new_chain_id) == 0) {
                            break;
                        }
                    }
                }
                
                // Add chain to target structure
                gemmi::Chain new_chain;
                new_chain.name = new_chain_id;
                for (const gemmi::Residue& res : source_chain.residues) {
                    new_chain.residues.push_back(res);
                }
                
                // Find or create model in target structure
                if (target_structure.models.empty()) {
                    target_structure.models.emplace_back();
                }
                target_structure.models[0].chains.push_back(new_chain);
                
                if (std::find(chain_ids.begin(), chain_ids.end(), new_chain_id) == chain_ids.end()) {
                    chain_ids.push_back(new_chain_id);
                }
            }
            status = true;
        }
    }

    return std::make_pair(status, chain_ids);
}

} // namespace coot