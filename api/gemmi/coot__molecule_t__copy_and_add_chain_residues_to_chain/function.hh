#pragma once
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace molecule_t {

// Copy residues from source_chain to target_chain, appending them
// Both chains must belong to different gemmi::Structure objects
void copy_and_add_chain_residues_to_chain_gemmi(
    const gemmi::Structure& source_structure,
    gemmi::Structure& target_structure,
    const std::string& source_chain_id,
    const std::string& target_chain_id);

} // namespace molecule_t
} // namespace coot