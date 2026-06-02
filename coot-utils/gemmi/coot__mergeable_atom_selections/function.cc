#include "function.hh"
#include <cmath>
#include <algorithm>

namespace coot {

// Implementation of match_container_t_gemmi::add
void match_container_t_gemmi::add(gemmi::Atom* at_1, gemmi::Residue* res_1,
                                   gemmi::Atom* at_2, gemmi::Residue* res_2) {
  for (auto& m : matches) {
    if (m.residue_1 == res_1 && m.residue_2 == res_2) {
      m.add(at_1, at_2);
      return;
    }
  }
  match_container_for_residues_t_gemmi new_match(res_1, res_2);
  new_match.add(at_1, at_2);
  matches.push_back(std::move(new_match));
}

// Implementation of match_container_t_gemmi::find_best_match
match_container_for_residues_t_gemmi
match_container_t_gemmi::find_best_match() const {
  match_container_for_residues_t_gemmi best;
  double best_dist_sq = 1e30;

  for (const auto& m : matches) {
    if (m.atom_pairs.size() >= 2) {
      double dist_sq = 0;
      for (const auto& p : m.atom_pairs) {
        double dx = p.first->pos.x - p.second->pos.x;
        double dy = p.first->pos.y - p.second->pos.y;
        double dz = p.first->pos.z - p.second->pos.z;
        dist_sq += dx*dx + dy*dy + dz*dz;
      }
      if (dist_sq < best_dist_sq) {
        best = m;
        best_dist_sq = dist_sq;
      }
    }
  }
  return best;
}

// Main function: gemmi port of mergeable_atom_selections
std::pair<bool, match_container_for_residues_t_gemmi>
mergeable_atom_selections_gemmi(
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel1,
    const std::vector<std::pair<gemmi::Atom*, gemmi::Residue*>>& sel2) {

  bool status = false;
  match_container_for_residues_t_gemmi m;

  double max_dist = 0.8;

  if (!sel1.empty() && !sel2.empty()) {
    match_container_t_gemmi match_set;

    for (const auto& e1 : sel1) {
      for (const auto& e2 : sel2) {
        gemmi::Atom* at_1 = e1.first;
        gemmi::Atom* at_2 = e2.first;

        // Skip self-contacts (identical atom pointers)
        if (at_1 == at_2) continue;

        std::string atom_name_1 = at_1->name;
        std::string atom_name_2 = at_2->name;

        if (atom_name_1 == atom_name_2) {
          if (atom_name_1 != " O  ") {
            double dx = at_1->pos.x - at_2->pos.x;
            double dy = at_1->pos.y - at_2->pos.y;
            double dz = at_1->pos.z - at_2->pos.z;
            double dist_sq = dx*dx + dy*dy + dz*dz;
            if (dist_sq <= max_dist * max_dist) {
              match_set.add(at_1, e1.second, at_2, e2.second);
            }
          }
        }
      }
    }

    match_container_for_residues_t_gemmi best_match = match_set.find_best_match();

    if (best_match.residue_1) {
      m = best_match;
      status = true;
    }
  }

  return std::pair<bool, match_container_for_residues_t_gemmi>(status, m);
}

}  // namespace coot