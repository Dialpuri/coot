#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__residue__residue__c6495e/gemmi/function.hh"

namespace coot {
namespace backrub {

coot::minimol::residue
make_residue_include_only_gemmi(const gemmi::Residue* orig_prev_residue,
                                const std::vector<std::string>& prev_res_atoms) {
    coot::minimol::residue r(orig_prev_residue, prev_res_atoms);
    return r;
}

} // namespace backrub
} // namespace coot