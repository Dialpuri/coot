#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {

// Ported from MMDB's deep_copy_this_residue_old_style.
// Takes a gemmi::CRA to provide parent context (chain/residue).
// Returns a newly allocated Residue (caller owns it), or nullptr if no atoms were copied.
gemmi::Residue *
deep_copy_this_residue_old_style_gemmi(const gemmi::CRA &cra,
                                       const std::string &altconf,
                                       int whole_residue_flag,
                                       bool embed_in_chain_flag);

} // namespace coot