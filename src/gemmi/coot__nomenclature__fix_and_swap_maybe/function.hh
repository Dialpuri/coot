#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {
class protein_geometry;
}

namespace coot {
namespace nomenclature {

std::vector<gemmi::CRA> fix_and_swap_maybe_gemmi(
    gemmi::Structure& st,
    const coot::protein_geometry* geom_p,
    bool apply_swaps);

} // namespace nomenclature
} // namespace coot