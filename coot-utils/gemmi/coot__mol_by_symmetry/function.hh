#pragma once

#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <vector>

namespace coot {

/// Apply a fractional symmetry operation (as gemmi::FTransform) to a deep copy
/// of the structure, converting through orthogonal space using the provided cell.
/// An optional lattice pre-shift (in fractional abc units) is applied before and
/// undone after the symmetry operation.
///
/// gemmi port of coot::mol_by_symmetry (originally used mmdb::Manager and
/// clipper::Cell / clipper::RTop_frac).
gemmi::Structure mol_by_symmetry_gemmi(
    const gemmi::Structure& mol,
    gemmi::UnitCell cell,
    gemmi::FTransform ftrans,
    std::vector<int> pre_shift_to_origin_abc);

} // namespace coot