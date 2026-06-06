#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace molecule_t {

    int cis_trans_conversion_gemmi(const std::string &atom_cid,
                                   gemmi::Structure &st,
                                   gemmi::Model *standard_residues_model);

} // molecule_t
} // coot