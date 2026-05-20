#pragma once
#include <string>
#include <vector>
#include <gemmi/cif.hpp>

namespace coot {
namespace protein_geometry {

void comp_tree_gemmi(const gemmi::cif::Loop& loop, int imol_enc);

} // namespace protein_geometry
} // namespace coot