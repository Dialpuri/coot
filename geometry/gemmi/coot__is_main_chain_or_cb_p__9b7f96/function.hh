#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {

/** Check if an atom name corresponds to a main-chain atom or CB.
 * gemmi port of coot::is_main_chain_or_cb_p (string overload).
 */
bool is_main_chain_or_cb_p_gemmi(const std::string &atom_name);

/** Check if a gemmi::Atom is a main-chain atom or CB.
 * gemmi port of coot::is_main_chain_or_cb_p (atom overload).
 */
bool is_main_chain_or_cb_p_gemmi(const gemmi::Atom &atom);

} // namespace coot
