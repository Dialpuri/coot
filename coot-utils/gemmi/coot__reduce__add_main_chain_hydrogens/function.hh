#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_main_chain_H/gemmi/function.hh"

namespace coot::reduce {

inline bool is_standard_amino_acid_name(const std::string& name) {
    static const std::vector<std::string> standard = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL"
    };
    for (const auto& s : standard) {
        if (name == s) return true;
    }
    return false;
}

void add_main_chain_hydrogens_gemmi(gemmi::Residue& residue, gemmi::Residue* residue_prev, bool is_gly);

} // namespace coot::reduce