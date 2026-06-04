#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <string>

namespace coot {
namespace util {

void transform_selection_gemmi(
    gemmi::Model& model,
    const std::string& chain_name,
    const clipper::RTop_orth& rtop);

} // namespace util
} // namespace coot