#pragma once
#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace molecule_t {

void replace_coords_gemmi(
    gemmi::Model& mol_model,
    const std::vector<gemmi::CRA>& selection,
    bool change_altconf_occs_flag,
    bool replace_coords_with_zero_occ_flag
);

} // namespace molecule_t
} // namespace coot