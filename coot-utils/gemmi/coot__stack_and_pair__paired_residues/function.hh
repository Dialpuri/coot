#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <algorithm>
#include <cmath>

#include <clipper/core/coords.h>

#include "geometry/protein-geometry.hh"
#include <coot/geometry/hb-types.hh>

namespace coot {

namespace stack_and_pair {

   // gemmi version of paired_residues_info_t
   class paired_residues_info_t {
   public:
      paired_residues_info_t(gemmi::Residue *r1, gemmi::Residue *r2,
                             const std::vector<std::pair<gemmi::Atom *, gemmi::Atom *> > &atom_pair_vec_in) :
         res_1(r1), res_2(r2), atom_pair_vec(atom_pair_vec_in) {}
      gemmi::Residue *res_1;
      gemmi::Residue *res_2;
      std::vector<std::pair<gemmi::Atom *, gemmi::Atom *> > atom_pair_vec;
   };

   // gemmi port of paired_residues
   std::vector<stack_and_pair::paired_residues_info_t>
   paired_residues_gemmi(gemmi::Model &model,
                         const std::vector<std::pair<bool, gemmi::Residue *>> &residues_vec,
                         bool residues_are_all_moving_flag,
                         const coot::protein_geometry &geom,
                         int imol_enc);

} // namespace stack_and_pair
} // namespace coot