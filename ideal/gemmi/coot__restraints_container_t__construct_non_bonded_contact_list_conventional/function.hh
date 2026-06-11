#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

namespace coot {

class restraints_container_t {
public:
  std::vector<gemmi::Atom*> gemmi_atoms;
  std::vector<std::vector<int>> bonded_atom_indices;
  std::vector<std::vector<int>> filtered_non_bonded_atom_indices;
  std::vector<std::vector<gemmi::Residue*>> active_residues_vec;
  int n_atoms = 0;
  bool have_oxt_flag = false;

  void filter_non_bonded_by_distance_gemmi(
    const std::vector<const gemmi::Atom*>& atoms,
    const std::vector<std::vector<int>>& non_bonded_atom_indices,
    double dist);

  void construct_non_bonded_contact_list_conventional_gemmi();
};

inline void coot::restraints_container_t::filter_non_bonded_by_distance_gemmi(
    const std::vector<const gemmi::Atom*>& atoms,
    const std::vector<std::vector<int>>& non_bonded_atom_indices,
    double dist)
{
    filtered_non_bonded_atom_indices.clear();
    filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());

    double dist_lim2 = dist * dist;

    for (unsigned int i = 0; i < non_bonded_atom_indices.size(); i++) {
        for (unsigned int j = 0; j < non_bonded_atom_indices[i].size(); j++) {
            const gemmi::Atom& atom_1 = *atoms[i];
            const gemmi::Atom& atom_2 = *atoms[non_bonded_atom_indices[i][j]];

            gemmi::Vec3 diff = atom_1.pos - atom_2.pos;
            double dist2 = diff.length_sq();

            if (dist2 < dist_lim2) {
                int i_at_ind = non_bonded_atom_indices[i][j];
                filtered_non_bonded_atom_indices[i].push_back(i_at_ind);
            }
        }
    }
}

inline void coot::restraints_container_t::construct_non_bonded_contact_list_conventional_gemmi() {

  // Helper to find index of an atom in gemmi_atoms
  auto find_atom_index = [this](const gemmi::Atom* atom) -> int {
    for (int i = 0; i < n_atoms; i++) {
      if (gemmi_atoms[i] == atom) return i;
    }
    return -1;
  };

  std::vector<std::vector<int>> non_bonded_atom_indices;
  non_bonded_atom_indices.resize(bonded_atom_indices.size());

  // Iterate over all active residues (flattened from active_residues_vec)
  std::vector<gemmi::Residue*> all_active_residues;
  for (auto& chain_residues : active_residues_vec) {
    for (auto* res : chain_residues) {
      all_active_residues.push_back(res);
    }
  }

  int nSelResidues_active = static_cast<int>(all_active_residues.size());

  for (int iar = 0; iar < nSelResidues_active; iar++) {
    gemmi::Residue* res = all_active_residues[iar];

    for (const gemmi::Atom& atom_ref : res->atoms) {
      const gemmi::Atom* atom_ptr = &atom_ref;

      int atom_index = find_atom_index(atom_ptr);
      if (atom_index < 0) continue;

      // Skip OXT if have_oxt_flag
      bool matched_oxt = false;
      if (have_oxt_flag) {
        if (atom_ptr->name == " OXT") {
          matched_oxt = true;
        }
      }
      if (matched_oxt) continue;

      // Compare with all atoms in all active residues
      for (int jar = 0; jar < nSelResidues_active; jar++) {
        gemmi::Residue* res_inner = all_active_residues[jar];

        for (const gemmi::Atom& atom_inner_ref : res_inner->atoms) {
          const gemmi::Atom* atom_inner_ptr = &atom_inner_ref;

          int atom_index_inner = find_atom_index(atom_inner_ptr);
          if (atom_index_inner < 0) continue;

          if (atom_index == atom_index_inner) {
            continue;
          }

          bool inner_matched_oxt = false;
          if (have_oxt_flag) {
            if (atom_inner_ptr->name == " OXT") {
              inner_matched_oxt = true;
            }
          }
          if (inner_matched_oxt) continue;

          // Check if bonded
          bool was_bonded_flag = false;
          for (int bonded_idx : bonded_atom_indices[atom_index]) {
            if (bonded_idx == atom_index_inner) {
              was_bonded_flag = true;
              break;
            }
          }

          if (!was_bonded_flag) {
            non_bonded_atom_indices[atom_index].push_back(atom_index_inner);
          }
        }
      }
    }
  }

  // Convert to const pointers for the filter call
  std::vector<const gemmi::Atom*> const_atoms;
  const_atoms.reserve(gemmi_atoms.size());
  for (auto* a : gemmi_atoms) {
    const_atoms.push_back(a);
  }

  filter_non_bonded_by_distance_gemmi(const_atoms, non_bonded_atom_indices, 8.0);
}

} // namespace coot
