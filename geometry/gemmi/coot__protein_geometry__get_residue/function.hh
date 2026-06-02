#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

// Free-function port of coot::protein_geometry::get_residue.
// Returns a newly-allocated gemmi::Residue (caller must delete).
gemmi::Residue* get_residue_gemmi(
    protein_geometry& geom,
    const std::string& comp_id,
    int imol_enc,
    bool idealised_flag,
    bool try_autoload_if_needed,
    float b_factor);

} // namespace coot