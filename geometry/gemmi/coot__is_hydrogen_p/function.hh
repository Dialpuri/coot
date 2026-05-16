#pragma once

#include <gemmi/model.hpp>

namespace coot {

// Port of coot::is_hydrogen_p — returns true if the atom's element is
// hydrogen (H) or deuterium (D).
inline bool is_hydrogen_p_gemmi(const gemmi::Atom *at) {
    return at->element.is_hydrogen();
}

}  // namespace coot