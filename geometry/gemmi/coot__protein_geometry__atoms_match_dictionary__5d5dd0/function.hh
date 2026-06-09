#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"

namespace coot {

std::pair<bool, std::vector<std::string>> atoms_match_dictionary_gemmi(
    const gemmi::Residue& residue,
    bool check_hydrogens_too_flag,
    bool apply_bond_distance_check) {

  std::string res_name(residue.name);

  coot::protein_geometry geom;
  geom.init_standard();

  std::pair<bool, dictionary_residue_restraints_t> restraints =
      geom.get_monomer_restraints(res_name, 0);

  if (!restraints.first) {
    return std::pair<bool, std::vector<std::string>>(false, std::vector<std::string>());
  }

  std::vector<std::string> atom_name_vec;
  bool status = true;

  for (const auto& atom : residue.atoms) {
    std::string residue_atom_name(coot::pad_atom_name(atom.name, atom.element.name()));

    bool found = false;

    // Skip hydrogen checking if flag is off (mirrors MMDB: ele == " H" || ele == " D")
    if (atom.element.is_hydrogen() && !check_hydrogens_too_flag) {
      found = true;
    }

    if (!found) {
      for (const auto& restraint_atom : restraints.second.atom_info) {
        if (restraint_atom.atom_id_4c == residue_atom_name) {
          found = true;
          atom_name_vec.push_back(residue_atom_name);
          break;
        }
      }
    }

    if (!found && residue_atom_name != " OXT") {
      if (std::find(atom_name_vec.begin(), atom_name_vec.end(), residue_atom_name) == atom_name_vec.end()) {
        atom_name_vec.push_back(residue_atom_name);
      }
      status = false;
    }
  }

  if (status && apply_bond_distance_check) {
    // Inline bond distance check
    auto should_include_atom = [&](const gemmi::Atom& a) {
      if (check_hydrogens_too_flag) return true;
      return !a.element.is_hydrogen();
    };

    for (size_t ibond = 0; ibond < restraints.second.bond_restraint.size(); ibond++) {
      const auto& bond = restraints.second.bond_restraint[ibond];
      std::string bond_atom_1 = bond.atom_id_1_4c();
      std::string bond_atom_2 = bond.atom_id_2_4c();

      for (size_t iat = 0; iat < residue.atoms.size(); iat++) {
        const gemmi::Atom& at_1 = residue.atoms[iat];
        std::string an1 = coot::pad_atom_name(at_1.name, at_1.element.name());
        if (an1 == bond_atom_1 && should_include_atom(at_1)) {
          for (size_t jat = 0; jat < residue.atoms.size(); jat++) {
            const gemmi::Atom& at_2 = residue.atoms[jat];
            std::string an2 = coot::pad_atom_name(at_2.name, at_2.element.name());
            if (an2 == bond_atom_2 && should_include_atom(at_2)) {
              if (at_1.altloc == at_2.altloc) {
                double dx = at_1.pos.x - at_2.pos.x;
                double dy = at_1.pos.y - at_2.pos.y;
                double dz = at_1.pos.z - at_2.pos.z;
                double d_sq = dx*dx + dy*dy + dz*dz;
                if (d_sq > 10.0) {
                  status = false;
                  break;
                }
              }
            }
          }
        }
        if (!status) break;
      }
      if (!status) break;
    }
  }

  return std::pair<bool, std::vector<std::string>>(status, atom_name_vec);
}

} // namespace coot