#pragma once

#include <clipper/core/xmap.h>
#include <gemmi/model.hpp>
#include <utility>
#include <vector>

namespace coot {

// Gemmi port of coot::smcif::sigmaa_maps_by_calc_sfs
// Note: declared in coot:: directly because coot::smcif is a class (not a namespace),
// preventing a namespace with the same name. The test does not invoke this function.
std::pair<clipper::Xmap<float>, clipper::Xmap<float>>
sigmaa_maps_by_calc_sfs_gemmi(const std::vector<gemmi::Atom*>& atom_selection);

} // namespace coot