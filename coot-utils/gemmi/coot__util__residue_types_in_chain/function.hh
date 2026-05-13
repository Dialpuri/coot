#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

// Original: coot::util::residue_types_in_chain(mmdb::Chain *chain_p)
// Returns unique residue types in a chain (sorted alphabetically due to set)
std::vector<std::string>
residue_types_in_chain_gemmi(const gemmi::Chain &chain);

} // namespace util
} // namespace coot