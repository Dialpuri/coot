#pragma once

#include <string>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_atom/gemmi/function.hh"

namespace coot {

float molecule_t_get_temperature_factor_of_atom_gemmi(const std::string &cid, gemmi::Structure &st);

} // namespace coot