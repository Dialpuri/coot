#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/core/coords.h>
#include <array>

namespace coot {

void transform_by_gemmi(gemmi::Structure& st, std::array<std::array<double, 4>, 4> mat);

} // namespace coot