#pragma once

#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include "ideal/extra-restraints.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

struct restraints_container_t {
   std::vector<std::pair<bool, void*>> residues_vec;
   bool from_residue_vector;
   std::vector<std::vector<int>> restraints_vec;
   std::map<atom_spec_t, int> atom_map;
   int udd_atom_index_handle;
   void* mol;
   std::set<int> fixed_atom_indices;

   bool fixed_check(int idx) const {
      return fixed_atom_indices.count(idx) > 0;
   }

   void add_user_defined_target_position_restraint(
       int type, int atom_index, const atom_spec_t& atom_spec,
       const clipper::Coord_orth& pos, double weight);
};

void add_extra_target_position_restraints_gemmi(restraints_container_t& rc,
                                                 const extra_restraints_t& extra_restraints);

} // namespace coot