#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>

namespace coot {

bool insert_coords_internal_gemmi(const gemmi::Model& src_model, gemmi::Model& target_model);

} // namespace coot