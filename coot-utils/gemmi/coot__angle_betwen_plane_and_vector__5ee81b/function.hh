#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <string>
#include <vector>
#include <utility>

namespace coot {

std::pair<bool, double> angle_betwen_plane_and_vector_gemmi(
    gemmi::Residue *residue_p,
    const std::vector<std::string> &ring_atom_names,
    const std::string &altconf_in,
    const clipper::Coord_orth &vector);

} // namespace coot