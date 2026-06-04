#pragma once

#include "ideal/simple-restraint.hh"
#include "geometry/residue-and-atom-specs.hh"

#include <string>
#include <vector>

namespace coot {

std::string format_gemmi(
    simple_restraint const& r,
    std::vector<atom_spec_t> const& atoms,
    double distortion);

} // namespace coot