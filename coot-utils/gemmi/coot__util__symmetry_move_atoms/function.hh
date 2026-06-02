#pragma once

#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <vector>
#include <utility>

namespace coot {
namespace util {

std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>>
symmetry_move_atoms_gemmi(
    const std::vector<clipper::Coord_orth>& protein_coords,
    const std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>>& water_atoms_in,
    clipper::Cell cell,
    clipper::Spacegroup spacegroup);

} // namespace util
} // namespace coot