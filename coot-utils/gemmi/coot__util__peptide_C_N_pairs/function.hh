#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

inline bool is_standard_amino_acid_name(const std::string& name) {
    static const std::string standard[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL"
    };
    for (const auto& s : standard) {
        if (s == name) return true;
    }
    return false;
}

std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>
peptide_C_N_pairs_gemmi(gemmi::Chain* chain_p);

} // namespace util
} // namespace coot