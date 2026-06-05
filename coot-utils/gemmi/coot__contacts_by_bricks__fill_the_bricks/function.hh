#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <set>

namespace coot {
namespace contacts_by_bricks {

void fill_the_bricks_gemmi(
    gemmi::Atom const* atoms[],
    int n_atoms,
    float brick_size,
    const float lower_left[3],
    const int range[3],
    std::vector<std::set<unsigned int>>& atoms_in_bricks);

} // namespace contacts_by_bricks
} // namespace coot