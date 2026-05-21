#pragma once

#include <string>
#include <clipper/clipper.h>

namespace coot {
namespace smcif {

clipper::Cell get_cell_for_data_gemmi(const std::string& file_name);

} // namespace smcif
} // namespace coot