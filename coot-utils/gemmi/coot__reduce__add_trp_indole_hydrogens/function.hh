#pragma once

#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_trp_indole_hydrogen/gemmi/function.hh"

namespace coot::reduce {

void add_trp_indole_hydrogens_gemmi(gemmi::Residue& residue, double bl_HN, double bl_arom) {
    add_trp_indole_hydrogen_gemmi(" HD1", " CG ", " CD1", " NE1", bl_HN, residue);
    add_trp_indole_hydrogen_gemmi(" HE1", " CD1", " NE1", " CE2", bl_arom, residue);
    add_trp_indole_hydrogen_gemmi(" HE3", " CD2", " CE3", " CZ3", bl_arom, residue);
    add_trp_indole_hydrogen_gemmi(" HZ3", " CE3", " CZ3", " CH2", bl_arom, residue);
    add_trp_indole_hydrogen_gemmi(" HH2", " CZ3", " CH2", " CZ2", bl_arom, residue);
    add_trp_indole_hydrogen_gemmi(" HZ2", " CH2", " CZ2", " CE2", bl_arom, residue);
}

} // namespace coot::reduce