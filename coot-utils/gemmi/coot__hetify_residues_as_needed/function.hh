#pragma once

#include <gemmi/model.hpp>

namespace coot {

/// Port of coot::hetify_residues_as_needed to gemmi.
/// Iterates over all chains and residues in the given model,
/// calling hetify_residue_atoms_as_needed_gemmi on each residue.
/// Always returns 0.
int hetify_residues_as_needed_gemmi(gemmi::Model* model);

} // namespace coot