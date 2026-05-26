#pragma once
#include <gemmi/pdb.hpp>
#include <string>
#include <vector>
#include <utility>

namespace coot {

/// Port of coot::molecule_t::try_add_by_consolidation to gemmi
/// Tries to consolidate chains from source_structure into target_structure
/// by matching single-residue-type chains. Returns (success, list_of_modified_chain_ids)
std::pair<bool, std::vector<std::string>>
molecule_t_try_add_by_consolidation_gemmi(const gemmi::Structure& source_structure,
                                          gemmi::Structure& target_structure);

} // namespace coot