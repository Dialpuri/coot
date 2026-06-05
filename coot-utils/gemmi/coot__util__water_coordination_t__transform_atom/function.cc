#include "function.hh"

namespace coot {
namespace util {

void water_coordination_t::transform_atom_gemmi(int i, int j) {
  contact_atoms_info_t& cai = atom_contacts[i];
  contact_atom_t& cat = cai[j];
  gemmi::Atom* at = cat.at;

  double x = at->pos.x;
  double y = at->pos.y;
  double z = at->pos.z;

  // Apply 4x4 homogeneous transformation matrix (row-major convention:
  // mat[0][3], mat[1][3], mat[2][3] hold translation)
  double new_x = cat.mat[0][0]*x + cat.mat[0][1]*y + cat.mat[0][2]*z + cat.mat[0][3];
  double new_y = cat.mat[1][0]*x + cat.mat[1][1]*y + cat.mat[1][2]*z + cat.mat[1][3];
  double new_z = cat.mat[2][0]*x + cat.mat[2][1]*y + cat.mat[2][2]*z + cat.mat[2][3];

  at->pos.x = static_cast<float>(new_x);
  at->pos.y = static_cast<float>(new_y);
  at->pos.z = static_cast<float>(new_z);
}

}  // namespace util
}  // namespace coot