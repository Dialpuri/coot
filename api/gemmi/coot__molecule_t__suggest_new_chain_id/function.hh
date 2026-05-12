#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <gemmi/model.hpp>

namespace coot {

// Suggest a new chain ID not already present in the structure.
// current_chain_id is the chain ID we're considering adding.
// Returns a chain ID that doesn't exist in the structure, or empty string if none found.
inline std::string molecule_t_suggest_new_chain_id_gemmi(const gemmi::Structure& st,
                                                         const std::string& current_chain_id) {
    std::string new_chain_id;
    
    // 20200110-PE no more specials at all
    // std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%^&@?/~|-+=(){}:;.,'");
    std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
    
    std::vector<std::string> existing;
    
    // Collect all existing chain IDs from the first model
    if (!st.models.empty()) {
        const gemmi::Model& model = st.models[0];
        for (const gemmi::Chain& chain : model.chains) {
            existing.push_back(chain.name);
        }
        
        // Check if current_chain_id exists
        bool found_it = false;
        for (const std::string& chid : existing) {
            if (chid == current_chain_id) {
                found_it = true;
                break;
            }
        }
        
        if (!found_it) {
            new_chain_id = current_chain_id; // all done!
        }
        
        // How about a multichar post-fix?
        if (new_chain_id.empty()) {
            if (current_chain_id.length() > 1) {
                std::string trial_chain_id = current_chain_id + "2";
                // mmdb chain id max length is 4 (mmdb_defs.h)
                if (trial_chain_id.length() < 4) {
                    bool found_it = false;
                    for (const std::string& chid : existing) {
                        if (chid == trial_chain_id) {
                            found_it = true;
                            break;
                        }
                    }
                    if (!found_it) {
                        new_chain_id = trial_chain_id; // all done!
                    }
                }
            }
        }
        
        if (new_chain_id.empty()) { // not set yet
            unsigned int l = r.length();
            std::vector<std::string> candidates(l);
            for (unsigned int i = 0; i < l; i++)
                candidates[i] = std::string(1, r[i]);
            
            for (const std::string& ext : existing) {
                candidates.erase(std::remove(candidates.begin(), candidates.end(), ext), candidates.end());
            }
            
            if (!candidates.empty())
                new_chain_id = candidates[0];
        }
    }
    
    return new_chain_id;
}

} // namespace coot