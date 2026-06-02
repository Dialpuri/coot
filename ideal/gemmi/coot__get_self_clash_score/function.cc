#include "function.hh"
#include <cmath>

namespace coot {

double get_self_clash_score_gemmi(
    const gemmi::Structure& st,
    gemmi::Model& model,
    const std::vector<gemmi::CRA>& atom_selection,
    const std::vector<torsion_atom_quad_gemmi>& quads)
{
  // Score is
  // sum of (d-bump_max)^2 for atom pairs i,j where j<i where d < bump_max

  const double bump_max = 2.8; // 20170615 try this (for less self bumping)
  double clash_score = 0.0;
  int n = static_cast<int>(atom_selection.size());

  // Brute-force pairwise comparison (mirrors original SeekContacts with id1 < id2)
  // This exactly mirrors the MMDB logic: iterate all pairs (i,j) with i<j,
  // filter same-residue, hydrogen, O5, torsion pairs, then accumulate (bump_max - dist)^2
  for (int i = 0; i < n; i++) {
    if (!atom_selection[i].atom) continue;
    for (int j = i + 1; j < n; j++) {
      if (!atom_selection[j].atom) continue;

      const gemmi::CRA& cra1 = atom_selection[i];
      const gemmi::CRA& cra2 = atom_selection[j];

      // Must be in different residues (mirrors: at_1->residue != at_2->residue)
      if (!cra1.residue || !cra2.residue) continue;
      if (cra1.residue == cra2.residue) continue;

      // Neither element should be hydrogen (mirrors: e1 != " H" && e2 != " H")
      if (cra1.atom->element.is_hydrogen() || cra2.atom->element.is_hydrogen())
        continue;

      // Ignore bumps to O5 — check only the second atom (mirrors: atom_name_2 != " O5 ")
      // In gemmi, atom.name is trimmed, so "O5" in gemmi corresponds to " O5 " in MMDB
      if (cra2.atom->name == "O5")
        continue;

      // Calculate squared distance
      double dx = cra1.atom->pos.x - cra2.atom->pos.x;
      double dy = cra1.atom->pos.y - cra2.atom->pos.y;
      double dz = cra1.atom->pos.z - cra2.atom->pos.z;
      double d_sqd = dx * dx + dy * dy + dz * dz;

      // Only consider pairs within bump_max distance (mirrors SeekContacts filter)
      if (d_sqd >= bump_max * bump_max)
        continue;

      // Are they in a torsion? (mirrors: both_in_a_torsion_p)
      bool in_a_tors = both_in_a_torsion_p_gemmi(cra1, cra2, quads);
      if (!in_a_tors) {
        double dist = std::sqrt(d_sqd);
        double delta = bump_max - dist;
        clash_score += delta * delta;
      }
    }
  }

  return clash_score;
}

} // namespace coot