#pragma once

#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

// gemmi port of coot::atom_overlaps_container_t::symmetry_contacts(float d)
// Returns pairs of (atom1*, atom2*) for contacts within distance d
// across symmetry operations (excluding the identity at lattice origin).
std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>
symmetry_contacts_gemmi(
    const residue_spec_t& res_central,
    gemmi::Structure* st,
    float d
);

} // namespace coot