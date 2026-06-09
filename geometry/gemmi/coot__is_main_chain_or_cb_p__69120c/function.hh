#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

/// Port of coot::is_main_chain_or_cb_p from MMDB to gemmi.
/// Returns true if the atom name corresponds to a main-chain atom (N, C, CA, O, OXT, H, HA)
/// or a CB atom. Accepts both MMDB-padded names (e.g. " CA ") and unpadded gemmi names (e.g. "CA").
bool is_main_chain_or_cb_p_gemmi(const std::string& atom_name);

/// Overload that takes a gemmi::Atom directly.
bool is_main_chain_or_cb_p_gemmi(const gemmi::Atom& atom);

} // namespace coot