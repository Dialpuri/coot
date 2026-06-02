#pragma once

#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

namespace dipole {

std::vector<std::pair<gemmi::Atom*, float>>
charged_atoms_gemmi(
   std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> dict_res_pairs);

} // namespace dipole

} // namespace coot