#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
namespace util {

// Return all residues in the chain.
// The clicked_residue parameter is accepted for API compatibility but ignored.
inline std::vector<gemmi::Residue*>
get_residues_in_fragment_gemmi(gemmi::Chain* chain,
                               const coot::residue_spec_t& clicked_residue) {
    std::vector<gemmi::Residue*> result;
    for (auto& residue : chain->residues) {
        result.push_back(&residue);
    }
    return result;
}

}
}