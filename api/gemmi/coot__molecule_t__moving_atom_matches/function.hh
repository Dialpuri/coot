#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {
namespace molecule_t {

// Compare a moving atom (as CRA) against a reference atom in a selection
// Returns true if name, ins_code, alt_loc, chain_id, and resno match
bool moving_atom_matches_gemmi(const std::vector<gemmi::CRA>& selection,
                               size_t ref_idx,
                               const gemmi::CRA& moving_cra);

} // namespace molecule_t
} // namespace coot