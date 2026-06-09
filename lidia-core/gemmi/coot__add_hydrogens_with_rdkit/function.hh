#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>
#include <vector>

namespace coot {
class dictionary_residue_restraints_t;
}

namespace coot {

std::pair<bool, std::string> add_hydrogens_with_rdkit_gemmi(
    gemmi::CRA cra,
    const coot::dictionary_residue_restraints_t &restraints);

} // namespace coot
