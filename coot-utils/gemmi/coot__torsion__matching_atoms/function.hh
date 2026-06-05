#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_spec_t__matches_spec/gemmi/function.hh"

namespace coot {
namespace torsion {

std::vector<gemmi::Atom*> matching_atoms_gemmi(
    gemmi::Chain* chain,
    gemmi::Residue* residue,
    const std::pair<int, coot::atom_spec_t>& atom_1,
    const std::pair<int, coot::atom_spec_t>& atom_2,
    const std::pair<int, coot::atom_spec_t>& atom_3,
    const std::pair<int, coot::atom_spec_t>& atom_4);

} // namespace torsion
} // namespace coot