#pragma once

#include <gemmi/cifdoc.hpp>
#include <string>

namespace coot {
namespace protein_geometry {

// gemmi port of chem_comp: reads a CIF Loop (e.g. _chem_comp) and registers entries
std::string chem_comp_gemmi(gemmi::cif::Loop& loop, int imol_enc);

} // namespace protein_geometry
} // namespace coot