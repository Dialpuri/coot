#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <utility>
#include <string>
#include <vector>

#include "coot/mini-mol/atom-quads.hh"

namespace coot {
namespace match_torsions {

enum torsion_type_enum {
  REFERENCE_TORSION = 0,
  MOVING_TORSION = 1
};

static const gemmi::Atom*
find_atom_by_name(const gemmi::Residue &res, const std::string &atom_name) {
  for (const auto &atom : res.atoms) {
    // gemmi atom names are space-padded to 4 chars; compare after trimming
    std::string trimmed = atom.name;
    trimmed.erase(0, trimmed.find_first_not_of(' '));
    trimmed.erase(trimmed.find_last_not_of(' ') + 1);
    if (trimmed == atom_name) {
      return &atom;
    }
  }
  return nullptr;
}

/* Find four atoms by name in a gemmi::Residue and compute their dihedral
 * angle.  Returns {true, angle_in_radians} on success or {false, 0} when
 * any of the four atoms is missing.
 */
inline std::pair<bool, double>
torsion_angle(const gemmi::Residue &res, const atom_name_quad &quad) {
  const gemmi::Atom *a = find_atom_by_name(res, quad.atom_name(0));
  const gemmi::Atom *b = find_atom_by_name(res, quad.atom_name(1));
  const gemmi::Atom *c = find_atom_by_name(res, quad.atom_name(2));
  const gemmi::Atom *d = find_atom_by_name(res, quad.atom_name(3));

  if (!a || !b || !c || !d) return {false, 0.0};

  double angle = gemmi::calculate_dihedral(a->pos, b->pos, c->pos, d->pos);
  return {true, angle};
}

/* gemmi port of match_torsions::get_torsion.
 *
 * torsion_type selects which residue to use (REFERENCE_TORSION or MOVING_TORSION).
 * quad provides the four atom names.
 */
inline std::pair<bool, double>
get_torsion_gemmi(int torsion_type,
                  const gemmi::Residue &res_ref,
                  const gemmi::Residue &res_moving,
                  const coot::atom_name_quad &quad) {
  switch (torsion_type) {
    case REFERENCE_TORSION:
      return torsion_angle(res_ref, quad);
    case MOVING_TORSION:
      return torsion_angle(res_moving, quad);
    default:
      return std::pair<bool, double>(false, 0.0);
  }
}

} // namespace match_torsions
} // namespace coot
