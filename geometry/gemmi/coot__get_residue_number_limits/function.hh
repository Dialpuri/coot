#pragma once

#include <map>
#include <string>
#include <gemmi/model.hpp>

namespace coot {

// Returns per-chain min/max residue sequence numbers.
// Takes a gemmi::Structure (no MMDB dependency).
std::map<std::string, std::pair<int, int> >
get_residue_number_limits_gemmi(const gemmi::Structure& st);

} // namespace coot