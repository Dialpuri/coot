#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

// Port of coot::restraints_container_t::add_link_angle to gemmi.
// Returns the number of link angle restraints found for the given link type
// between two residues.
int add_link_angle_gemmi(std::string link_type,
                         gemmi::CRA first, gemmi::CRA second,
                         short int is_fixed_first,
                         short int is_fixed_second,
                         const coot::protein_geometry &geom);

} // namespace coot