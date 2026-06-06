#pragma once
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__switch_his_protonation__ee6360/gemmi/function.hh"

namespace coot {
namespace reduce {

// 2-argument wrapper: sets bl_arom = 1.082 (copied from add_hydrogen_atoms)
// then delegates to the 3-argument _gemmi port.
void switch_his_protonation_gemmi(gemmi::Residue* residue_p,
                                   gemmi::Atom*   current_H_atom) {
    double bl_arom = 1.082; // copied from add_hydrogen_atoms()
    switch_his_protonation_gemmi(residue_p, current_H_atom, bl_arom);
}

} // namespace reduce
} // namespace coot