#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
namespace molecule_t {

// Gemmi port of coot::molecule_t::full_atom_spec_to_atom_index
// Looks up an atom by specification (chain, residue number, insertion code,
// atom name, alt conf) in the first model of the structure.
// Returns the 0-based atom index within the model, or -1 if not found.
int full_atom_spec_to_atom_index_gemmi(
    const coot::atom_spec_t &atom_spec,
    gemmi::Structure &st);

} // namespace molecule_t
} // namespace coot