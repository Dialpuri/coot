#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

bool move_atom_gemmi(const std::string &atom_name_in, gemmi::Residue *res_p, const gemmi::Vec3 &new_pos);

} // namespace coot