#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "analysis/mogul-interface.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

coot::dictionary_residue_restraints_t
make_restraints_gemmi(
    const gemmi::Residue* residue_p,
    const std::string& comp_id,
    const std::vector<coot::mogul_item>& items,
    int imol,
    const coot::protein_geometry& geom
);

} // namespace coot