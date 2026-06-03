#pragma once

#include <gemmi/model.hpp>

namespace coot {
namespace contacts_by_bricks {

// gemmi port of set_lower_left_and_range
// Computes the lower-left corner and brick range for a set of atoms.
//
// Parameters:
//   atoms       - array of pointers to gemmi atoms
//   n_atoms     - number of atoms
//   brick_size  - size of each brick (original default: 20.0)
//   lower_left  - output: minimum coordinates minus 6.0
//   range       - output: number of bricks along each axis (+1)
void set_lower_left_and_range_gemmi(
    gemmi::Atom const* atoms[],
    int n_atoms,
    float brick_size,
    float lower_left[3],
    int range[3]);

} // namespace contacts_by_bricks
} // namespace coot