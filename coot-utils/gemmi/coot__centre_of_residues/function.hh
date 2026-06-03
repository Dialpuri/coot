#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>

namespace coot {

std::pair<bool, gemmi::Vec3>
centre_of_residues_gemmi(const std::vector<gemmi::Residue*> &residues);

} // namespace coot