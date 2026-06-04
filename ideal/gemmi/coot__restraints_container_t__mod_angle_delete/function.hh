#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>

// Restraint type constant (matches coot::simple_restraint.hh)
constexpr int ANGLE_RESTRAINT_TYPE = 2;

// ----------------------------------------------------------------
// chem_mod_angle_gemmi — self-contained gemmi version of chem_mod_angle
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

// ----------------------------------------------------------------
// restraints_container_t_gemmi — self-contained gemmi version
// ----------------------------------------------------------------
struct restraints_container_t_gemmi {
  // Each restraint stored as [type, idx1, idx2, idx3]
  std::vector<std::vector<int>> restraints_vec;
  // bonded_atom_indices[i] = set of atoms bonded to atom i
  std::vector<std::set<int>> bonded_atom_indices;
  // Maps global atom index → residue pointer (for matching during delete)
  std::vector<gemmi::Residue*> atom_residues;
  // Maps global atom index → atom name (trimmed, for matching during delete)
  std::vector<std::string> atom_names;

  int size() const { return static_cast<int>(restraints_vec.size()); }

  void mod_angle_delete_gemmi(const chem_mod_angle_gemmi &mod_angle,
                              gemmi::Residue *residue_p);
};

// ----------------------------------------------------------------
// mod_angle_delete_gemmi implementation (preserves original bug:
// name_2 == mod_angle.atom_id_3 instead of name_3)
// ----------------------------------------------------------------
inline void restraints_container_t_gemmi::mod_angle_delete_gemmi(
    const chem_mod_angle_gemmi &mod_angle,
    gemmi::Residue *residue_p)
{
  std::vector<std::vector<int>>::iterator it;
  for (it = restraints_vec.begin(); it != restraints_vec.end(); ) {
    if (it->size() >= 4 && (*it)[0] == ANGLE_RESTRAINT_TYPE) {
      int idx1 = (*it)[1];
      int idx2 = (*it)[2];
      int idx3 = (*it)[3];
      if (static_cast<size_t>(idx1) < atom_residues.size() &&
          atom_residues[idx1] == residue_p) {
        if (static_cast<size_t>(idx2) < atom_residues.size() &&
            atom_residues[idx2] == residue_p) {
          std::string name_1 = (static_cast<size_t>(idx1) < atom_names.size())
                                  ? atom_names[idx1] : "";
          std::string name_2 = (static_cast<size_t>(idx2) < atom_names.size())
                                  ? atom_names[idx2] : "";
          std::string name_3 = (static_cast<size_t>(idx3) < atom_names.size())
                                  ? atom_names[idx3] : "";
          if (name_1 == mod_angle.atom_id_1) {
            if (name_2 == mod_angle.atom_id_2) {
              if (name_2 == mod_angle.atom_id_3) { // BUG: should be name_3
                it = restraints_vec.erase(it);
                continue;
              }
            }
          }
        }
      }
    }
    ++it;
  }
}