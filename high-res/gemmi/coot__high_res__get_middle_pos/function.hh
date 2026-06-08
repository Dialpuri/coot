#pragma once

#include <clipper/clipper.h>
#include <gemmi/model.hpp>
#include <utility>

namespace coot {
namespace high_res {

std::pair<clipper::Coord_orth, gemmi::Structure*>
get_middle_pos_gemmi(const gemmi::Structure& st);

} // namespace high_res
} // namespace coot
