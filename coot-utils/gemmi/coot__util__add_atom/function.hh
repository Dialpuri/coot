#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>
#include <string>

namespace coot {
namespace util {

bool add_atom_gemmi(
    gemmi::Residue* res,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const std::string& atom_name_3,
    const std::string& alt_conf,
    double length,
    double angle,
    double torsion,
    const std::string& new_atom_name,
    const std::string& new_atom_ele,
    float new_atom_occ,
    float new_atom_b_factor
);

} // namespace util
} // namespace coot