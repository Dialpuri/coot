#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <unordered_set>

namespace coot {

// Gemmi-compatible torsion quad: stores gemmi::CRA pointers instead of mmdb::Atom*
struct torsion_atom_quad_gemmi {
  gemmi::CRA atom_1{nullptr, nullptr, nullptr};
  gemmi::CRA atom_2{nullptr, nullptr, nullptr};
  gemmi::CRA atom_3{nullptr, nullptr, nullptr};
  gemmi::CRA atom_4{nullptr, nullptr, nullptr};
  double angle = 0;
  double deviation = 0;
  bool use = false;
  int atom_1_bond_order = 1;
  int atom_4_bond_order = 1;

  torsion_atom_quad_gemmi() = default;

  bool has_4_atoms() const {
    return atom_1.atom != nullptr && atom_2.atom != nullptr &&
           atom_3.atom != nullptr && atom_4.atom != nullptr;
  }

  bool is_valid() const {
    return has_4_atoms();
  }
};

// Gemmi version of both_in_a_torsion_p
// Returns true if both atoms are in any of the torsion quads
inline bool both_in_a_torsion_p_gemmi(
    const gemmi::CRA& a1,
    const gemmi::CRA& a2,
    const std::vector<torsion_atom_quad_gemmi>& quads)
{
  for (const auto& q : quads) {
    if (!q.is_valid()) continue;
    bool atom1 = (q.atom_1.atom == a1.atom) ||
                 (q.atom_2.atom == a1.atom) ||
                 (q.atom_3.atom == a1.atom) ||
                 (q.atom_4.atom == a1.atom);
    bool atom2 = (q.atom_1.atom == a2.atom) ||
                 (q.atom_2.atom == a2.atom) ||
                 (q.atom_3.atom == a2.atom) ||
                 (q.atom_4.atom == a2.atom);
    if (atom1 && atom2)
      return true;
  }
  return false;
}

// Port of coot::get_self_clash_score to gemmi API
// Uses gemmi::NeighborSearch for contact finding
double get_self_clash_score_gemmi(
    const gemmi::Structure& st,
    gemmi::Model& model,
    const std::vector<gemmi::CRA>& atom_selection,
    const std::vector<torsion_atom_quad_gemmi>& quads);

} // namespace coot