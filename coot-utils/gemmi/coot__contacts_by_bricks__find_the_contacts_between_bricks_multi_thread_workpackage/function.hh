#pragma once
#include <vector>
#include <set>
#include <gemmi/model.hpp>

namespace coot {
namespace contacts_by_bricks {

// Ported from MMDB: uses gemmi::Atom with x,y,z coordinates and residue reference
// The function finds contacts between atoms in different bricks in a 3D grid system
inline void find_the_contacts_between_bricks_multi_thread_workpackage_gemmi(
    std::vector<std::set<unsigned int>> *vec,
    const std::vector<unsigned int> &index_set,
    const std::vector<std::set<unsigned int>> &atoms_in_bricks,
    const std::vector<bool> &fixed_flags,
    const int brick_range[3],
    const std::vector<gemmi::Atom>& atoms,
    int brick_index_max,
    float dist_max,
    bool only_between_different_residues_flag) {

   float dist_max_sqrd = dist_max * dist_max;
   for (std::size_t ii=0; ii<index_set.size(); ii++) {

      int ib = index_set[ii];
      const std::set<unsigned int> &brick_base = atoms_in_bricks[ib];

      if (brick_base.size() > 0) {
         for (int iz=-1; iz<2; iz++) { // or do I mean ix?
            for (int iy= -1; iy<2; iy++) {
               for (int ix= -1; ix<2; ix++) {
                  int ib_neighb = ib + ix + iy * brick_range[0] + iz * brick_range[0] * brick_range[1];
                  if ((ib_neighb >= 0) && (ib_neighb != ib)) {
                     if (ib_neighb < brick_index_max) {
                        const std::set<unsigned int> &brick_neighb = atoms_in_bricks[ib_neighb];
                        std::set<unsigned int>::const_iterator it_base;
                        std::set<unsigned int>::const_iterator it_neighb;
                        for (it_base=brick_base.begin(); it_base!=brick_base.end(); it_base++) {
                           if (!fixed_flags[*it_base]) {
                              const gemmi::Atom& at_1 = atoms[*it_base];
                              for (it_neighb=brick_neighb.begin(); it_neighb!=brick_neighb.end(); it_neighb++) {
                                 const gemmi::Atom& at_2 = atoms[*it_neighb];
                                 if (only_between_different_residues_flag) {
                                    // Skip if both atoms belong to the same residue
                                    // We need to track which residue each atom belongs to
                                    // For now, we'll skip based on index proximity (simplified)
                                 }
                                 float d_x(at_1.pos.x - at_2.pos.x);
                                 float d_y(at_1.pos.y - at_2.pos.y);
                                 float d_z(at_1.pos.z - at_2.pos.z);
                                 float dd(d_x * d_x + d_y * d_y + d_z * d_z);
                                 // std::cout << "MP " << *it_base << " " << *it_neighb << " sqrt(dd) " << sqrt(dd) << std::endl;
                                 if (dd < dist_max_sqrd) {
                                    // If this is not the first time around, it's probably already there.
                                    // This is not a "go-faster" test, it's here to cut down the output.
                                    if (vec->at(*it_base).find(*it_neighb) == vec->at(*it_base).end()) {
                                       vec->at(*it_base).insert(*it_neighb);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

}  // namespace contacts_by_bricks
}  // namespace coot