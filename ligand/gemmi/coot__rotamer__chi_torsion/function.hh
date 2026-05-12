#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <clipper/core/clipper_util.h>

namespace coot {
namespace rotamer {

double chi_torsion_gemmi(const std::vector<int> &chi_angle_atoms_indices,
                         const std::vector<gemmi::Atom>& residue_atoms);

} // namespace rotamer
} // namespace coot