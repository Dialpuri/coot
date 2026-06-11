#pragma once

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include <gemmi/model.hpp>
#include <string>

namespace coot::util {

clipper::Xmap<float> make_gaussian_atom_map_for_mask_gemmi(
    const clipper::Xmap<float> &map_ref,
    gemmi::Model &mol,
    const std::string &cid,
    float sigma,
    float box_radius);

} // namespace coot::util
