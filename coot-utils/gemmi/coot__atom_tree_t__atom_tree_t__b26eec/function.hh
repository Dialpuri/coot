#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <gemmi/model.hpp>
#include <stdexcept>

#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"
#include "geometry/protein-geometry.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_name_map/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__atom_tree_t__fill_atom_vertex_vec_using_contacts/gemmi/function.hh"

namespace coot {

class atom_tree_t {
public:
  int made_from_minimol_residue_flag;
  int n_selected_atoms;
  int* atom_selection;

  std::map<std::string, map_index_t> name_map;
  std::vector<std::vector<int>> contact_indices;
  std::vector<gemmi::Atom> residue_atoms;
  std::vector<atom_vertex> atom_vertex_vec;
  std::vector<std::pair<int, int>> bonds;

  // Gemmi port of the original constructor
  void atom_tree_t_gemmi(
      const coot::dictionary_residue_restraints_t &rest,
      const gemmi::Residue &res,
      const std::string &altconf) {

    made_from_minimol_residue_flag = 1;
    n_selected_atoms = 0;
    atom_selection = nullptr;

    // Copy atoms from the gemmi residue, filtering by altconf
    residue_atoms.clear();
    char alt = (altconf.empty()) ? '\0' : altconf[0];
    for (const auto &atom : res.atoms) {
      if (atom.altloc == alt) {
        residue_atoms.push_back(atom);
      }
    }

    // Build name map using existing gemmi port
    name_map = fill_name_map_gemmi(res, altconf);

    // Build contact indices from restraints
    contact_indices.clear();
    size_t n_atoms = residue_atoms.size();
    contact_indices.resize(n_atoms);
    for (size_t i = 0; i < n_atoms; i++) {
      contact_indices[i].clear();
    }

    // Look up bonded atom pairs from restraints (bond_restraint holds dict_bond_restraint_t)
    for (const auto &bond : rest.bond_restraint) {
      std::string name_a = bond.atom_id_1_4c();
      std::string name_b = bond.atom_id_2_4c();

      int idx_a = -1;
      int idx_b = -1;

      if (name_map.count(name_a)) {
        idx_a = name_map.at(name_a).index();
      }
      if (name_map.count(name_b)) {
        idx_b = name_map.at(name_b).index();
      }

      if (idx_a >= 0 && idx_b >= 0) {
        contact_indices[idx_a].push_back(idx_b);
        contact_indices[idx_b].push_back(idx_a);
      }
    }

    // Build atom vertex vec using existing gemmi port
    atom_vertex_vec.clear();
    bool ok = fill_atom_vertex_vec_using_contacts_gemmi(
        contact_indices, 0, residue_atoms, atom_vertex_vec);

    if (!ok) {
      throw std::runtime_error("Failed to construct atom tree");
    }
  }
};

} // namespace coot