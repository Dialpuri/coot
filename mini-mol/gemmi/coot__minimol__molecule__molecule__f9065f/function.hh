#pragma once

#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "mini-mol/mini-mol.hh"

namespace coot {
namespace minimol {

/// gemmi port of coot::minimol::molecule default constructor.
/// Returns a default-constructed molecule (have_cell=0, have_spacegroup=0, no fragments).
inline molecule molecule_gemmi() {
    return molecule{};
}

} // namespace minimol
} // namespace coot
