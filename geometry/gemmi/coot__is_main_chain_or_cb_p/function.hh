#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {
  // Port of coot::is_main_chain_or_cb_p from MMDB to gemmi
  // Returns true if atom is a main chain atom or CB atom
  bool is_main_chain_or_cb_p_gemmi(const gemmi::Atom& atom);

  // Overload that takes atom name string
  bool is_main_chain_or_cb_p_gemmi(const std::string& atom_name);
}