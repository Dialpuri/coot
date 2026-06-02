#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "coot/mini-mol/mini-mol.hh"

namespace coot {
namespace residue_by_phi_psi {

clipper::Coord_orth best_fit_phi_psi_attaching_oxygen_position_update_gemmi(
    const minimol::molecule &mm, const gemmi::Residue *residue_p);

} // namespace residue_by_phi_psi
} // namespace coot