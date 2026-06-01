#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {

inline bool atom_name_is(const gemmi::Atom& atom, const std::string& target) {
  // gemmi atom.name can be padded (" C  ") or unpadded ("C") depending on source.
  // Strip spaces for comparison.
  std::string stripped = atom.name;
  stripped.erase(stripped.find_last_not_of(" \t\r\n") + 1);
  stripped.erase(0, stripped.find_first_not_of(" \t\r\n"));
  return stripped == target;
}

inline std::pair<bool, bool>
peptide_C_and_N_are_close_p_gemmi(const gemmi::Residue* r1, const gemmi::Residue* r2) {
    float dist_crit = 2.8f;

    const gemmi::Atom* at_c_1 = nullptr;
    const gemmi::Atom* at_n_1 = nullptr;
    const gemmi::Atom* at_c_2 = nullptr;
    const gemmi::Atom* at_n_2 = nullptr;

    for (const auto& atom : r1->atoms) {
        if (atom_name_is(atom, "C")) { at_c_1 = &atom; }
        if (atom_name_is(atom, "N")) { at_n_1 = &atom; }
    }

    for (const auto& atom : r2->atoms) {
        if (atom_name_is(atom, "C")) { at_c_2 = &atom; }
        if (atom_name_is(atom, "N")) { at_n_2 = &atom; }
    }

    if (at_c_1 && at_n_2) {
        double d = (at_c_1->pos - at_n_2->pos).length();
        if (d < dist_crit)
            return {true, false};
    }

    if (at_n_1 && at_c_2) {
        double d = (at_n_1->pos - at_c_2->pos).length();
        if (d < dist_crit)
            return {true, true};
    }

    return {false, false};
}

} // namespace coot