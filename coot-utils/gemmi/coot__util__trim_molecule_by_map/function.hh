#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/xmap.h"

namespace coot {
namespace util {

// Mirror the enum values from coot-trim.hh (avoid including it since it references mmdb)
constexpr int TRIM_BY_MAP_DELETE  = 0;
constexpr int TRIM_BY_MAP_ZERO_OCC = 1;

int trim_molecule_by_map_gemmi(gemmi::Model &model,
                               const clipper::Xmap<float> &xmap,
                               float map_level,
                               short int remove_or_zero_occ_flag,
                               short int waters_only_flag);

} // namespace util
} // namespace coot
