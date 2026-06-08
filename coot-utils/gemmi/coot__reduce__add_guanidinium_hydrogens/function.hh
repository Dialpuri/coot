#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>
#include <iostream>

#include <clipper/core/coords.h>

// Include actual generated dependency headers
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bond_length_angle_torsion/gemmi/function.hh"

namespace coot {
namespace reduce {

void add_guanidinium_hydrogens_gemmi(gemmi::Residue& residue);

} // namespace reduce
} // namespace coot
