#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {
namespace util {

class contact_atom_t {
 public:
  double dist;
  double mat[4][4];
  gemmi::Atom *at;

  contact_atom_t(gemmi::Atom *contactor, gemmi::Atom *central_atom)
    : dist(0), mat{{}}, at(contactor)
  {
    gemmi::Vec3 c1(contactor->pos.x, contactor->pos.y, contactor->pos.z);
    gemmi::Vec3 c2(central_atom->pos.x, central_atom->pos.y, central_atom->pos.z);
    dist = (c1 - c2).length();
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        mat[i][j] = (i == j) ? 1.0 : 0.0;
  }

  contact_atom_t(gemmi::Atom *contactor, gemmi::Atom *central_atom, const double mat_in[4][4])
    : dist(0), at(contactor)
  {
    gemmi::Vec3 c1(contactor->pos.x, contactor->pos.y, contactor->pos.z);
    gemmi::Vec3 c2(central_atom->pos.x, central_atom->pos.y, central_atom->pos.z);
    dist = (c1 - c2).length();
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        mat[i][j] = mat_in[i][j];
  }
};

class contact_atoms_info_t {
 public:
  gemmi::Atom *central_atom;
  std::vector<contact_atom_t> contactors;

  contact_atoms_info_t(gemmi::Atom *central, contact_atom_t contactor)
    : central_atom(central)
  {
    contactors.push_back(std::move(contactor));
  }

  contact_atom_t& operator[](int index) {
    return contactors[index];
  }
};

class water_coordination_t {
 public:
  std::vector<contact_atoms_info_t> atom_contacts;

  void transform_atom_gemmi(int i, int j);
};

}  // namespace util
}  // namespace coot