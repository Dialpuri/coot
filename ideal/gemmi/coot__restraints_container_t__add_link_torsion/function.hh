#pragma once

#include <gemmi/model.hpp>
#include <string>

#include "geometry/protein-geometry.hh"

namespace coot {

// Port of coot::restraints_container_t::add_link_torsion to gemmi.
// Iterates through link torsion restraints matching the given link_type,
// finds matching atoms in the two residues, skips mainchain torsions
// (N-N, CA-CA, C-C), and counts them.
// Returns 0 (matching the original function's behavior of never incrementing n_torsions).
int add_link_torsion_gemmi(std::string link_type,
                           gemmi::CRA first,
                           gemmi::CRA second,
                           short int is_fixed_first,
                           short int is_fixed_second,
                           const coot::protein_geometry &geom);

} // namespace coot