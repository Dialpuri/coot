#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>

namespace coot {

// Port of have_intermediate_residue_by_seqnum to gemmi.
// Takes gemmi::CRA (Chain*, Residue*, Atom*) for the two residues and
// a vector of gemmi::CRA entries to search through for an intermediate
// residue in the same chain.
bool have_intermediate_residue_by_seqnum_gemmi(
    gemmi::CRA first,
    gemmi::CRA second,
    const std::vector<std::pair<bool, gemmi::CRA>>& residues_vec);

}  // namespace coot