#include "function.hh"
#include <algorithm>
#include <gemmi/polyheur.hpp>

namespace coot {
namespace molecule_t {

void copy_and_add_chain_residues_to_chain_gemmi(
    const gemmi::Structure& source_structure,
    gemmi::Structure& target_structure,
    const std::string& source_chain_id,
    const std::string& target_chain_id) {
    
    // Find the source chain in source_structure
    const gemmi::Chain* source_chain = nullptr;
    for (const gemmi::Model& model : source_structure.models) {
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
                               [&source_chain_id](const gemmi::Chain& c) { return c.name == source_chain_id; });
        if (it != model.chains.end()) {
            source_chain = &*it;
            break;
        }
    }
    if (!source_chain) return;  // source chain not found
    
    // Find the target chain in target_structure
    gemmi::Chain* target_chain = nullptr;
    for (gemmi::Model& model : target_structure.models) {
        auto it = std::find_if(model.chains.begin(), model.chains.end(),
                               [&target_chain_id](const gemmi::Chain& c) { return c.name == target_chain_id; });
        if (it != model.chains.end()) {
            target_chain = &*it;
            break;
        }
    }
    if (!target_chain) return;  // target chain not found
    
    // Copy each residue from source chain to target chain
    for (const gemmi::Residue& res : source_chain->residues) {
        // Create a copy of the residue
        gemmi::Residue new_res = res;
        target_chain->residues.push_back(std::move(new_res));
    }
}

} // namespace molecule_t
} // namespace coot