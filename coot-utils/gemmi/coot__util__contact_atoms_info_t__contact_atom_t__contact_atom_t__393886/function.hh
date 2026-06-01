#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

namespace coot {
namespace util {

class contact_atoms_info_t {

 public:

  class contact_atom_t {

   public:
    double dist;
    double mat[4][4];
    gemmi::Atom const *at;

    // Default: identity matrix
    contact_atom_t(gemmi::Atom const *contactor,
                   gemmi::Atom const *central_atom)
      : dist(0), mat{{}}, at(contactor)
    {
      gemmi::Vec3 c1(contactor->pos.x, contactor->pos.y, contactor->pos.z);
      gemmi::Vec3 c2(central_atom->pos.x, central_atom->pos.y, central_atom->pos.z);
      dist = (c1 - c2).length();
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          mat[i][j] = (i == j) ? 1.0 : 0.0;
    }

    // With custom matrix
    contact_atom_t(gemmi::Atom const *contactor,
                   gemmi::Atom const *central_atom,
                   const double mat_in[4][4])
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

};

}  // namespace util
}  // namespace coot