#pragma once

#include "mini-mol/mini-mol.hh"
#include "clipper/core/coords.h"
#include <gemmi/model.hpp>

// Gemmi port of coot::high_res::fill_globular_protein
// Takes the globular molecule to fill as a parameter (by reference),
// the source molecule with atom data, the target position,
// and the gemmi::Structure for periodic boundary closest-approach.
void coot_high_res_fill_globular_protein_gemmi(
    coot::minimol::molecule &globular_mol,
    const coot::minimol::molecule &mol,
    const clipper::Coord_orth &target_pos_in,
    const gemmi::Structure &structure);