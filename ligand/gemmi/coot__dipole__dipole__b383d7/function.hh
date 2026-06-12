#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__dipole__init/gemmi/function.hh"

namespace coot {
namespace dipole {

dipole_state_t dipole_gemmi(
    const coot::dictionary_residue_restraints_t &rest,
    gemmi::CRA cra);

} // namespace dipole
} // namespace coot
