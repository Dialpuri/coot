#pragma once
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__is_nucleotide/gemmi/function.hh"

namespace coot { namespace stack_and_pair {

// Return true if any residue in the input vector is a nucleotide (RNA or DNA)
// This mirrors the original MMDB function logic that checked atoms and their residues
inline bool
contains_nucleic_acid_gemmi(const std::vector<gemmi::Residue>& residues) {
    for (const auto& residue : residues) {
        if (coot::util::is_nucleotide_gemmi(residue)) {
            return true;
        }
    }
    return false;
}

}} // namespace coot::stack_and_pair