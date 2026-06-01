#pragma once

#include <gemmi/model.hpp>

namespace coot {

/// Check that residues in a chain have progressive sequence numbers
/// (each residue's seqid.num.value >= previous + 1).
/// Returns 1 if progressive, 0 otherwise.
short int progressive_residues_in_chain_check_gemmi(const gemmi::Chain& chain);

}  // namespace coot