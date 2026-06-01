#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot { namespace molecule_t {

std::vector<std::string> get_residue_types_without_dictionaries_gemmi(
    const gemmi::Model &model,
    const coot::protein_geometry &geom,
    int imol_no);

}} // namespace coot::molecule_t