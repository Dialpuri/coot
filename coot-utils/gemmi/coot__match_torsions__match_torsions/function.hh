#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"

namespace coot {

// Data structure holding the gemmi-compatible state of match_torsions
// (mirrors the private member variables of the original MMDB class)
struct match_torsions_data_t {
   gemmi::Residue *res_moving;
   gemmi::Residue *res_ref;
   dictionary_residue_restraints_t moving_residue_restraints;

   match_torsions_data_t()
      : res_moving(nullptr),
        res_ref(nullptr) {}
};

// Ported version of coot::match_torsions::match_torsions constructor
// Stores gemmi residue pointers and restraints.
inline match_torsions_data_t match_torsions_gemmi(
   gemmi::Residue *res_moving_in,
   gemmi::Residue *res_ref_in,
   const dictionary_residue_restraints_t &rest
) {
   match_torsions_data_t data;
   data.res_moving = res_moving_in;
   data.res_ref = res_ref_in;
   data.moving_residue_restraints = rest;
   return data;
}

} // namespace coot