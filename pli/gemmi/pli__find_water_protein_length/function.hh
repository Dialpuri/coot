#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <gemmi/model.hpp>
#include "clipper/core/coords.h"

namespace pli {

// Find the minimum distance between atoms in a ligand residue and
// atoms in protein residues (excluding water molecules).
// Returns the distance (not squared).
double find_water_protein_length_gemmi(const gemmi::Residue* ligand_residue,
                                       const gemmi::Model* model);

} // namespace pli