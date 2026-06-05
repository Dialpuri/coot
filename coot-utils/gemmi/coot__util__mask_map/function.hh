#pragma once

#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <vector>

namespace coot { namespace util {

clipper::Xmap<float>
mask_map_gemmi(const clipper::Xmap<float> &xmap_in,
               const std::vector<const gemmi::Residue *> &neighb_residues);

}} // namespace coot::util