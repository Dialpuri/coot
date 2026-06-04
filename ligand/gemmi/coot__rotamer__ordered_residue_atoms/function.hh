#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace rotamer {

/// Return atoms of a residue in the canonical rotamer order.
/// For residues with exchangeable pairs (VAL CG1/CG2, etc.) this ensures
/// atom_order[0] always precedes atom_order[1].
/// For PRO the full backbone+sidechain order is enforced.
/// For all other residues atoms are returned in PDB order.
std::vector<gemmi::Atom*> ordered_residue_atoms_gemmi(const gemmi::Residue& residue);

} // namespace rotamer
} // namespace coot