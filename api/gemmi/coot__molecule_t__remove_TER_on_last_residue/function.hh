#pragma once
#include <gemmi/model.hpp>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-dev/coot/api/gemmi/coot__molecule_t__remove_TER_internal/function.hh"

namespace coot { namespace molecule_t {

void remove_TER_on_last_residue_gemmi(gemmi::Chain& chain) {
    if (!chain.residues.empty()) {
        gemmi::Residue& last_residue = chain.residues.back();
        remove_TER_internal_gemmi(last_residue);
    }
}

}}