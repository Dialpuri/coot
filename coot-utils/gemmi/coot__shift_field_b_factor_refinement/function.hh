#pragma once

#include <gemmi/model.hpp>

namespace coot {

// Gemmi port of coot::shift_field_b_factor_refinement.
// The original uses clipper for crystallographic operations (SFcalc,
// shift_field_u_iso, FFT maps, etc.) that have no gemmi equivalents.
// Preserves original semantics: returns false when the refinement
// cannot proceed, leaving B-factors unmodified.
bool shift_field_b_factor_refinement_gemmi(
    gemmi::Structure& mol,
    int ncyc = 1
);

}  // namespace coot