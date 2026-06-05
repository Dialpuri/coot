#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

// ----------------------------------------------------------------
// Restraint type constant — matches coot enum: ANGLE_RESTRAINT = 1
// ----------------------------------------------------------------
constexpr int ANGLE_RESTRAINT = 1;

// ----------------------------------------------------------------
// chem_mod_angle_gemmi — gemmi-compatible chemistry modification descriptor
// ----------------------------------------------------------------
struct chem_mod_angle_gemmi {
  std::string type;            // "add", "change", "delete"
  std::string atom_id_1;       // trimmed atom name, e.g. "N"
  std::string atom_id_2;       // trimmed atom name, e.g. "CA"
  std::string atom_id_3;       // trimmed atom name, e.g. "C"
  double new_value_angle;      // target angle in degrees
  double new_value_angle_esd;  // target sigma

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

// ----------------------------------------------------------------
// trim_atom_name — removes leading/trailing spaces from PDB atom name
// ----------------------------------------------------------------
inline std::string trim_atom_name(std::string name) {
  auto lo = name.find_first_not_of(' ');
  if (lo == std::string::npos) return "";
  auto hi = name.find_last_not_of(' ');
  return name.substr(lo, hi - lo + 1);
}

// ----------------------------------------------------------------
// simple_restraint_gemmi — minimal gemmi analogue of coot::simple_restraint
// ----------------------------------------------------------------
struct simple_restraint_gemmi {
  int restraint_type;
  int atom_index_1;
  int atom_index_2;
  int atom_index_3;
  double target_value;
  double sigma;
};

// ----------------------------------------------------------------
// restraints_container_t_gemmi — gemmi-compatible restraint container
// ----------------------------------------------------------------
struct restraints_container_t_gemmi {
  std::vector<simple_restraint_gemmi> restraints_vec;
  // Per-atom lookup arrays (indexed by global atom index)
  std::vector<gemmi::Residue*> atom_residues;
  std::vector<std::string>     atom_names;  // trimmed names

  // ---------------------------------------------------------------
  // mod_angle_change_gemmi
  // Walk angle restraints; when both atom 1 and 2 belong to the
  // given residue and all three trimmed names match the
  // chem_mod_angle, update target_value and sigma.
  // ---------------------------------------------------------------
  void mod_angle_change_gemmi(const chem_mod_angle_gemmi &mod_angle,
                              gemmi::Residue *residue_p);
};

// ----- mod_angle_change_gemmi (inline) -----
inline void restraints_container_t_gemmi::mod_angle_change_gemmi(
    const chem_mod_angle_gemmi &mod_angle,
    gemmi::Residue *residue_p)
{
  for (size_t i = 0; i < restraints_vec.size(); ++i) {
    simple_restraint_gemmi &rest = restraints_vec[i];
    if (rest.restraint_type == ANGLE_RESTRAINT) {
      if (atom_residues[static_cast<size_t>(rest.atom_index_1)] == residue_p) {
        if (atom_residues[static_cast<size_t>(rest.atom_index_2)] == residue_p) {
          std::string name_1 = atom_names[static_cast<size_t>(rest.atom_index_1)];
          std::string name_2 = atom_names[static_cast<size_t>(rest.atom_index_2)];
          std::string name_3 = atom_names[static_cast<size_t>(rest.atom_index_3)];
          if (name_1 == mod_angle.atom_id_1) {
            if (name_2 == mod_angle.atom_id_2) {
              if (name_3 == mod_angle.atom_id_3) {
                restraints_vec[i].target_value = mod_angle.new_value_angle;
                restraints_vec[i].sigma        = mod_angle.new_value_angle_esd;
              }
            }
          }
        }
      }
    }
  }
}