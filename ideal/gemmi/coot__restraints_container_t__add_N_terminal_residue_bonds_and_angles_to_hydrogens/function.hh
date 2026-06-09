#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <coot/ideal/simple-restraint.hh>

namespace coot {

// Port of coot::restraints_container_t::add_N_terminal_residue_bonds_and_angles_to_hydrogens
// Standalone free function that analyzes N-terminal residue hydrogen atoms
// (H1/H2/H3) and counts the bond+angle restraints that would be added.

inline coot::restraints_container_t::restraint_counts_t
add_N_terminal_residue_bonds_and_angles_to_hydrogens_gemmi(
    gemmi::Residue *residue_p
) {
  coot::restraints_container_t::restraint_counts_t result_counts;

  // Build a local atom-index map: gemmi::Atom* -> sequential int
  std::map<gemmi::Atom*, int> atom_index_map;
  int next_index = 0;
  for (gemmi::Atom &atom : residue_p->atoms) {
    atom_index_map[&atom] = next_index++;
  }

  // Inline N and CA index lookups (same logic as get_N_index / get_CA_index)
  // gemmi normalizes atom names (strips spaces), so match "N" and "CA"
  int N_index = -2;
  for (int i = 0; i < static_cast<int>(residue_p->atoms.size()); i++) {
    if (residue_p->atoms[i].name == "N") {
      N_index = i;
      break;
    }
  }

  int CA_index = -2;
  for (int i = 0; i < static_cast<int>(residue_p->atoms.size()); i++) {
    if (residue_p->atoms[i].name == "CA") {
      CA_index = i;
      break;
    }
  }

  // Maps from altloc string -> local atom index for H1/H2/H3
  std::map<std::string, int> h1s;
  std::map<std::string, int> h2s;
  std::map<std::string, int> h3s;

  // Track bonded pairs for counting
  std::set<std::pair<int,int>> bonded_pairs;

  for (gemmi::Atom &atom : residue_p->atoms) {
    std::string atom_name = atom.name;
    // gemmi normalizes atom names (strips spaces), so match "H1", "H2", "H3"
    if (atom_name == "H1" || atom_name == "H2" || atom_name == "H3") {
      if (N_index >= 0 && CA_index >= 0) {
        int atom_index_1 = atom_index_map[&atom];
        int atom_index_2 = atom_index_map[&(residue_p->atoms[static_cast<size_t>(N_index)])];
        int atom_index_3 = atom_index_map[&(residue_p->atoms[static_cast<size_t>(CA_index)])];

        // Count: 1 bond restraint (H-N) + 1 angle restraint (H-N-CA)
        result_counts.n_bond_restraints++;
        result_counts.n_angle_restraints++;

        // Track bonded indices
        bonded_pairs.insert({atom_index_1, atom_index_2});
        bonded_pairs.insert({atom_index_2, atom_index_1});
        bonded_pairs.insert({atom_index_1, atom_index_3});
        bonded_pairs.insert({atom_index_3, atom_index_1});
      }

      // Store atoms for inter-hydrogen angle restraints (altloc as key)
      std::string altloc_str(1, atom.altloc);
      if (atom_name == "H1") {
        h1s[altloc_str] = atom_index_map[&atom];
      }
      if (atom_name == "H2") {
        h2s[altloc_str] = atom_index_map[&atom];
      }
      if (atom_name == "H3") {
        h3s[altloc_str] = atom_index_map[&atom];
      }
    }
  }

  // Inter-hydrogen angle restraints
  if (N_index >= 0) {
    int N_ud_index = atom_index_map[&(residue_p->atoms[static_cast<size_t>(N_index)])];

    for (const auto &h1_pair : h1s) {
      const std::string &key_alt_conf = h1_pair.first;
      auto it_2 = h2s.find(key_alt_conf);
      auto it_3 = h3s.find(key_alt_conf);

      if (it_2 != h2s.end()) {
        // H1-N-H2 angle restraint
        result_counts.n_angle_restraints++;
        bonded_pairs.insert({h1_pair.second, it_2->second});
        bonded_pairs.insert({it_2->second, h1_pair.second});
      }

      if (it_3 != h3s.end()) {
        // H1-N-H3 angle restraint
        result_counts.n_angle_restraints++;
        bonded_pairs.insert({h1_pair.second, it_3->second});
        bonded_pairs.insert({it_3->second, h1_pair.second});
      }

      if (it_2 != h2s.end() && it_3 != h3s.end()) {
        // H2-N-H3 angle restraint
        result_counts.n_angle_restraints++;
        bonded_pairs.insert({it_2->second, it_3->second});
        bonded_pairs.insert({it_3->second, it_2->second});
      }
    }
  }

  return result_counts;
}

} // namespace coot
