#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include <vector>
#include <string>

namespace coot {
namespace util {

std::vector<std::vector<int>>
get_contact_indices_from_restraints_gemmi(gemmi::Residue *residue,
                                          const coot::dictionary_residue_restraints_t &restraints,
                                          bool regular_residue_flag,
                                          bool add_reverse_contacts);

} // namespace util
} // namespace coot