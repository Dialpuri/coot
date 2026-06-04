#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

namespace coot {

// ----------------------------------------------------------------
// Standalone types for the gemmi port
// ----------------------------------------------------------------

struct chem_mod_angle_gemmi {
  std::string type;
  std::string atom_id_1;
  std::string atom_id_2;
  std::string atom_id_3;
  double new_value_angle;
  double new_value_angle_esd;

  chem_mod_angle_gemmi(std::string t,
                       std::string a1, std::string a2, std::string a3,
                       double ang, double esd)
    : type(std::move(t)),
      atom_id_1(std::move(a1)),
      atom_id_2(std::move(a2)),
      atom_id_3(std::move(a3)),
      new_value_angle(ang),
      new_value_angle_esd(esd) {}
};

// Restraint type constant (matches coot::simple_restraint.hh)
constexpr int ANGLE_RESTRAINT_TYPE = 2;

// ----------------------------------------------------------------
// restraints_container_t_gemmi — self-contained gemmi version
// Mirrors the original coot::restraints_container_t internals.
// ----------------------------------------------------------------

struct restraints_container_t_gemmi {
  // Each restraint stored as [type, idx1, idx2, idx3]
  std::vector<std::vector<int>> restraints_vec;
  // bonded_atom_indices[i] = set of atoms bonded to atom i
  std::vector<std::set<int>> bonded_atom_indices;

  int size() const { return static_cast<int>(restraints_vec.size()); }

  // ----------------------------------------------------------------
  // make_fixed_flags — returns 3 booleans (one per atom)
  // ----------------------------------------------------------------
  std::vector<bool> make_fixed_flags(int, int, int) const {
    return std::vector<bool>(3, false);
  }

  // ----------------------------------------------------------------
  // is_hydrogen — checks if an atom is hydrogen
  // ----------------------------------------------------------------
  static bool is_hydrogen(const gemmi::Atom &atom) {
    return atom.element.is_hydrogen();
  }

  // ----------------------------------------------------------------
  // add — push a new restraint only if it doesn't already exist
  // (matches original coot behaviour where duplicate restraints are skipped)
  // ----------------------------------------------------------------
  bool add(int rest_type,
           int idx1, int idx2, int idx3,
           const std::vector<bool> &,
           double, double,
           bool)
  {
    // Check for duplicate: skip if identical restraint already present
    for (const auto &r : restraints_vec) {
      if (r.size() >= 4 &&
          r[0] == rest_type &&
          r[1] == idx1 &&
          r[2] == idx2 &&
          r[3] == idx3) {
        return false;  // duplicate — no-op, matches original
      }
    }

    std::vector<int> restraint;
    restraint.push_back(rest_type);
    restraint.push_back(idx1);
    restraint.push_back(idx2);
    restraint.push_back(idx3);
    restraints_vec.push_back(restraint);

    int max_idx = std::max({idx1, idx2, idx3});
    if (static_cast<size_t>(max_idx + 1) > bonded_atom_indices.size())
      bonded_atom_indices.resize(max_idx + 1);

    return true;
  }

  // ----------------------------------------------------------------
  // mod_angle_add_gemmi — the ported function
  // ----------------------------------------------------------------
  void mod_angle_add_gemmi(const chem_mod_angle_gemmi &mod_angle,
                           gemmi::Residue &residue);
};

// ----------------------------------------------------------------
// Helper: trim PDB atom name (removes padding spaces)
// ----------------------------------------------------------------
inline std::string trim_atom_name(std::string name) {
  auto lo = name.find_first_not_of(' ');
  if (lo == std::string::npos) return "";
  auto hi = name.find_last_not_of(' ');
  return name.substr(lo, hi - lo + 1);
}

// ----------------------------------------------------------------
// mod_angle_add_gemmi implementation
// ----------------------------------------------------------------
inline void restraints_container_t_gemmi::mod_angle_add_gemmi(
    const chem_mod_angle_gemmi &mod_angle,
    gemmi::Residue &residue)
{
  size_t n_atoms = residue.atoms.size();

  for (size_t iat_1 = 0; iat_1 < n_atoms; ++iat_1) {
    std::string name_1 = trim_atom_name(residue.atoms[iat_1].name);
    if (name_1 == mod_angle.atom_id_1) {
      for (size_t iat_2 = 0; iat_2 < n_atoms; ++iat_2) {
        std::string name_2 = trim_atom_name(residue.atoms[iat_2].name);
        if (name_2 == mod_angle.atom_id_2) {
          for (size_t iat_3 = 0; iat_3 < n_atoms; ++iat_3) {
            std::string name_3 = trim_atom_name(residue.atoms[iat_3].name);
            if (name_3 == mod_angle.atom_id_3) {
              // Check alt-loc compatibility
              // In gemmi, ' ' means "no alt loc" (equivalent to MMDB "")
              char alt_1 = residue.atoms[iat_1].altloc;
              char alt_2 = residue.atoms[iat_2].altloc;
              char alt_3 = residue.atoms[iat_3].altloc;

              if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
                  ((alt_1 == ' ') && (alt_2 == alt_3)) ||
                  ((alt_2 == ' ') && (alt_1 == alt_3)) ||
                  ((alt_3 == ' ') && (alt_1 == alt_2))) {

                int index_1 = static_cast<int>(iat_1);
                int index_2 = static_cast<int>(iat_2);
                int index_3 = static_cast<int>(iat_3);

                std::vector<bool> fixed_flags =
                    make_fixed_flags(index_1, index_2, index_3);

                bool is_single_Hydrogen_atom_angle_restraint = false;
                unsigned int nH = 0;
                if (is_hydrogen(residue.atoms[iat_1])) ++nH;
                if (is_hydrogen(residue.atoms[iat_3])) ++nH;
                if (nH == 1)
                  is_single_Hydrogen_atom_angle_restraint = true;

                add(ANGLE_RESTRAINT_TYPE, index_1, index_2, index_3,
                    fixed_flags,
                    mod_angle.new_value_angle,
                    mod_angle.new_value_angle_esd,
                    is_single_Hydrogen_atom_angle_restraint);

                // Ensure bonded_atom_indices is sized properly
                int max_idx = std::max({index_1, index_2, index_3});
                if (static_cast<size_t>(max_idx + 1) > bonded_atom_indices.size())
                  bonded_atom_indices.resize(max_idx + 1);

                bonded_atom_indices[index_1].insert(index_3);
                bonded_atom_indices[index_3].insert(index_1);
              }
            }
          }
        }
      }
    }
  }
}

} // namespace coot