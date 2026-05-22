#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <vector>

namespace coot {
namespace db_strands {

// Helper function to check if an atom is a main chain atom or CB
inline bool is_main_chain_or_cb_p_gemmi(const std::string& atom_name) {
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

inline bool is_main_chain_or_cb_p_gemmi(const gemmi::Atom& atom) {
  return is_main_chain_or_cb_p_gemmi(atom.name);
}

// Trim a gemmi Structure to main chain atoms only (including CB if present)
// This is the gemmi version of the original MMDB function
inline void trim_to_mainchain_gemmi(gemmi::Structure& st) {
    // Process each model in the structure
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                // Iterate backwards to safely remove atoms
                for (int iat = static_cast<int>(residue.atoms.size()) - 1; iat >= 0; --iat) {
                    gemmi::Atom& at = residue.atoms[iat];
                    std::string ele = at.element.name();
                    if (!is_main_chain_or_cb_p_gemmi(at)
                        || ele == "H" || ele == "D") {
                        residue.atoms.erase(residue.atoms.begin() + iat);
                    }
                }
            }
        }
    }
    // PDBCleanup and FinishStructEdit are not needed in gemmi
    // The structure is already in a valid state after atom removal
}

} // namespace db_strands
} // namespace coot