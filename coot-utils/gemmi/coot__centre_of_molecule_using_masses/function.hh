#pragma once

#include <string>
#include <map>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot {

std::pair<bool, clipper::Coord_orth>
centre_of_molecule_using_masses_gemmi(const gemmi::Structure& st);

} // namespace coot