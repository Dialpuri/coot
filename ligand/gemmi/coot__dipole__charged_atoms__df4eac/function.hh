#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>
#include "geometry/protein-geometry.hh"

namespace coot::dipole {

std::vector<std::pair<gemmi::Atom *, float>>
charged_atoms_gemmi(
    gemmi::Residue *residue_p,
    const coot::dictionary_residue_restraints_t &rest);

} // namespace coot::dipole