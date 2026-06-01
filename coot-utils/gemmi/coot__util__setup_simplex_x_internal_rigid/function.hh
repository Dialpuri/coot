#pragma once

#include <gsl/gsl_vector.h>
#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace util {

// gemmi port of setup_simplex_x_internal_rigid
// Original was a no-op: "All rots and trans angles are set to 0."
inline void setup_simplex_x_internal_rigid_gemmi(gsl_vector *x,
                                                 std::vector<gemmi::Atom*> atom_selection,
                                                 int n_selected_atoms) {
  // unnecessary, I think.  All rots and trans angles are set to 0.
}

}
}