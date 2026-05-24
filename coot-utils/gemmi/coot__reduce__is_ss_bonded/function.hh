#pragma once
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__hack_ss_bond_test/gemmi/function.hh"

namespace coot {
namespace reduce {

inline bool is_ss_bonded_gemmi(const gemmi::Residue* residue_p, const gemmi::Model& model) {
    bool status = false;
    if (residue_p) {
        if (residue_p->name == "CYS") {
            // Check SS bonds here
            bool i = hack_ss_bond_test_gemmi(residue_p, model);
            if (i) status = true;
        }
    }
    return status;
}

} // namespace reduce
} // namespace coot