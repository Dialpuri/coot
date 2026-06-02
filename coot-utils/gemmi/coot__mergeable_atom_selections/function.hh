#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <vector>
#include <string>

namespace coot {

// gemmi-equivalent of match_container_for_residues_t
// Groups matching atom pairs by their two residues
struct match_container_for_residues_t_gemmi {
  gemmi::Residue* residue_1;
  gemmi::Residue* residue_2;
  std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> atom_pairs;

  match_container_for_residues_t_gemmi() : residue_1(nullptr), residue_2(nullptr) {}
  match_container_for_residues_t_gemmi(gemmi::Residue* r1, gemmi::Residue* r2)
      : residue_1(r1), residue_2(r2) {}

  void add(gemmi::Atom* at_1, gemmi::Atom* at_2) {
    atom_pairs.push_back(std::make_pair(at_1, at_2));
  }
};

// gemmi-equivalent of match_container_t
// Collects match_container_for_residues_t for all residue-pair groups
struct match_container_t_gemmi {
  std::vector<match_container_for_residues_t_gemmi> matches;

  void add(gemmi::Atom* at_1, gemmi::Residue* res_1,
           gemmi::Atom* at_2, gemmi::Residue* res_2);

  match_container_for_residues_t_gemmi find_best_match() const;
};

// gemmi port of coot::mergeable_atom_selections
// Takes two selections of (atom, residue) pairs and determines if they
// have enough matching atom pairs at the same residue to merge.
// Each selection element is std::pair<gemmi::Atom*, gemmi::Residue*>.
std::pair<bool, match_container_for_residues_t_gemmi>
mergeable_atom_selections_gemmi(
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2);

}  // namespace coot