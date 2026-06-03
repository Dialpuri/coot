#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

// Transforms selected atoms by:
//   1. Translating by (-us,-vs,-ws) in fractional coordinates (mol_to_origin)
//   2. Applying symmetry operation symm_no
//   3. Adding lattice translation (symm_x,symm_y,symm_z) in fractional
//
// Returns transformed atom copies.
std::vector<gemmi::Atom> trans_sel_gemmi(
    gemmi::Structure* st,
    std::vector<gemmi::Atom> atoms,
    int symm_no,
    int symm_x, int symm_y, int symm_z,
    int us, int vs, int ws);

} // namespace coot