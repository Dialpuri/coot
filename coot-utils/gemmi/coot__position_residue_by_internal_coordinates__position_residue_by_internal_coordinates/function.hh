#pragma once

#include <gemmi/model.hpp>
#include "mini-mol/atom-quads.hh"

namespace coot {

// Forward declaration — see function.cc for include of actual port
gemmi::Atom *position_residue_by_internal_coordinates_get_atom_gemmi(
    gemmi::Residue *res_1, gemmi::Residue *res_2,
    const atom_name_quad &quad, int atom_index);

class position_residue_by_internal_coordinates_gemmi {
public:
  position_residue_by_internal_coordinates_gemmi(
    gemmi::Residue *residue_ref,
    gemmi::Residue *residue_moving,
    const atom_name_quad &quad,
    const double &bond_length,
    const double &bond_angle,
    const double &bond_torsion);

  bool move_moving_residue();

private:
  gemmi::Residue *residue_ref_;
  gemmi::Residue *residue_moving_;
  atom_name_quad quad_;
  double bond_length_;
  double bond_angle_;
  double bond_torsion_;
  gemmi::Atom *at_0_;
  gemmi::Atom *at_1_;
  gemmi::Atom *at_2_;
  gemmi::Atom *at_3_;
};

} // namespace coot