#pragma once
#include <vector>
#include <clipper/clipper.h>
#include <gemmi/model.hpp>

namespace coot {

std::vector<gemmi::Residue*>
residues_near_position(const clipper::Coord_orth &pt,
                       const gemmi::Structure &st,
                       double radius);

} // namespace coot