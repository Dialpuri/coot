#include "function.hh"
#include <cmath>
#include <gemmi/math.hpp>

namespace {

// Rewrite of coot::util::rotate_around_vector using gemmi::Vec3
// Rotate point `p` around axis `dir` by `angle` radians, pivoting at `base`
gemmi::Vec3 rotate_around_vector_gemmi(const gemmi::Vec3 &dir,
                                        const gemmi::Vec3 &p,
                                        const gemmi::Vec3 &base,
                                        double angle) {
  // Compute the direction vector and its length
  gemmi::Vec3 u = dir;
  double len = u.length();
  if (len < 1e-10) {
    // Degenerate axis - return original point
    return p;
  }
  u = u / len; // normalize

  // Vector from base to point
  gemmi::Vec3 v = p - base;

  // Parallel component (along u)
  double dot = u.dot(v);
  gemmi::Vec3 parallel = u * dot;

  // Perpendicular component
  gemmi::Vec3 perp = v - parallel;

  // Cross product gives direction of rotation
  gemmi::Vec3 cross = u.cross(perp);

  // Rotated position
  gemmi::Vec3 rotated = parallel + perp * std::cos(angle) + cross * std::sin(angle);

  return base + rotated;
}

} // anonymous namespace

void coot::atom_tree_t::rotate_internal_gemmi(
    std::vector<gemmi::Atom> &atoms,
    std::vector<coot::map_index_t> moving_atom_indices,
    const gemmi::Vec3 &dir,
    const gemmi::Vec3 &base_atom_pos,
    double angle) {

  for (unsigned int im = 0; im < moving_atom_indices.size(); im++) {
    int idx = moving_atom_indices[im].index();
    if (idx >= 0 && idx < static_cast<int>(atoms.size())) {
      gemmi::Atom &at = atoms[idx];
      gemmi::Vec3 po(at.pos.x, at.pos.y, at.pos.z);
      gemmi::Vec3 pt = rotate_around_vector_gemmi(dir, po, base_atom_pos, angle);
      at.pos.x = pt.x;
      at.pos.y = pt.y;
      at.pos.z = pt.z;
    }
  }
}