#pragma once

#include "clipper/core/cell.h"
#include "clipper/core/coords.h"
#include <string>

namespace coot {
namespace smcif {

clipper::Resolution get_resolution_gemmi(const clipper::Cell &cell,
                                         const std::string &file_name);

} // namespace smcif
} // namespace coot