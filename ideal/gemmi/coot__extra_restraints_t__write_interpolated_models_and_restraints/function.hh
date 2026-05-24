#pragma once

#include <string>
#include <map>
#include <gemmi/model.hpp>
#include "ideal/extra-restraints.hh"

namespace coot {

// Ported from MMDB version to use gemmi
void write_interpolated_models_and_restraints_gemmi(
    const extra_restraints_t &this_obj,
    const extra_restraints_t &final,
    const gemmi::Structure &st_1,
    const gemmi::Structure &st_2,
    unsigned int n_path_points,
    std::string file_name_stub);

} // namespace coot