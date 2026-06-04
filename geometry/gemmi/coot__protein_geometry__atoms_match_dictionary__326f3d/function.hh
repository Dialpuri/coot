#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__protein_geometry__atoms_match_dictionary__c725f6/gemmi/function.hh"

namespace coot {

std::pair<bool, std::vector<std::pair<gemmi::Residue*, std::vector<std::string>>>>
atoms_match_dictionary_gemmi(
    const std::vector<gemmi::Residue*>& residues,
    bool check_hydrogens_too_flag,
    bool apply_bond_distance_check) {

    bool status = true;
    std::vector<std::pair<gemmi::Residue*, std::vector<std::string>>> bad_residues;

    for (gemmi::Residue* res : residues) {
        auto r = coot::atoms_match_dictionary_gemmi(*res, check_hydrogens_too_flag, apply_bond_distance_check);
        if (!r.first) {
            bad_residues.push_back({res, r.second});
            status = false;
        }
    }

    return {status, bad_residues};
}

} // namespace coot