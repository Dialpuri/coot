#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__hole__assign_vdw_radii/gemmi/function.hh"

namespace coot::hole {

struct hole_result {
    gemmi::Model* model;
    clipper::Coord_orth from_pt;
    clipper::Coord_orth to_pt;
    double colour_map_multiplier;
    double colour_map_offset;
    int radius_handle;
};

inline hole_result hole_gemmi(
    gemmi::Model& model,
    clipper::Coord_orth from_pt_in,
    clipper::Coord_orth to_pt_in,
    const coot::protein_geometry& geom) {

    assign_vdw_radii_gemmi(model, geom);

    hole_result result;
    result.model = &model;
    result.from_pt = from_pt_in;
    result.to_pt = to_pt_in;
    result.colour_map_multiplier = 1.0;
    result.colour_map_offset = 0.0;
    result.radius_handle = 16777217;  // no UDR in gemmi; match original value
    return result;
}

} // namespace coot::hole