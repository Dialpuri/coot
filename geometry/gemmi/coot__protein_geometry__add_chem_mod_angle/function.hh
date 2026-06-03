#pragma once

#include <gemmi/cifdoc.hpp>
#include <string>

namespace coot {
class protein_geometry;
}

namespace coot {

void add_chem_mod_angle_gemmi(protein_geometry& geom, const gemmi::cif::Loop& loop);

} // namespace coot