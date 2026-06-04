#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__rotamer__probability_of_this_rotamer/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__deep_copy_residue/gemmi/function.hh"

namespace coot {
namespace rotamer_helper {

int optimize_rotamer_by_atom_names_gemmi(gemmi::CRA cra, bool apply_swap_if_found);

}  // namespace rotamer_helper
}  // namespace coot