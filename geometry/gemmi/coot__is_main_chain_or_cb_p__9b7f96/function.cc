#include "function.hh"
#include <string>

namespace coot {

bool is_main_chain_or_cb_p_gemmi(const std::string& atom_name) {
  // gemmi atom names are stored without padding, but may be padded in some contexts
  if (atom_name == " N  " ||
      atom_name == " C  " ||
      atom_name == " H  " ||
      atom_name == " CA " ||
      atom_name == " OXT" ||
      atom_name == " CB " ||
      atom_name == " HA " ||
      atom_name == " O  " ||
      atom_name == "N" ||
      atom_name == "C" ||
      atom_name == "H" ||
      atom_name == "CA" ||
      atom_name == "OXT" ||
      atom_name == "CB" ||
      atom_name == "HA" ||
      atom_name == "O") {
    return true;
  }
  return false;
}

bool is_main_chain_or_cb_p_gemmi(const gemmi::Atom& atom) {
  return is_main_chain_or_cb_p_gemmi(atom.name);
}

} // namespace coot
