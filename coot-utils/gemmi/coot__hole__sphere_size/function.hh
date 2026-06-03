#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <string>
#include <vector>

namespace coot {
    class protein_geometry;  // forward-declare coot type (defined in "geometry/protein-geometry.hh")
}

namespace coot::hole {

/// Gemmi port of coot::hole::sphere_size
///
/// Computes the radius of the largest sphere that can be centered at `pt`
/// without overlapping any atoms in `model` (restricted to the specified
/// chain names, or all chains if `chain_names` is empty).
///
/// For each atom, computes distance(pt, atom_pos) - vdw_radius.
/// Returns the minimum such value across all atoms, or -1 if there are no atoms.
///
/// Replaces the MMDB version that used selection handles and UDData.
double sphere_size_gemmi(
    const gemmi::Model& model,
    const std::vector<std::string>& chain_names,
    const gemmi::Vec3& pt,
    const coot::protein_geometry& geom
);

} // namespace coot::hole