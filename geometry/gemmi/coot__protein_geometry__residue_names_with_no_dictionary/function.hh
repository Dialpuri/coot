#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"

namespace coot {

std::vector<std::string>
residue_names_with_no_dictionary_gemmi(
    const protein_geometry& geom,
    const gemmi::Structure* st,
    int imol_no);

} // namespace coot