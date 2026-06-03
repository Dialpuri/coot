#pragma once

#define LIBCOOTAPI_BUILD
#include <rdkit/GraphMol/RDKitBase.h>
#include "pyrogen/restraints.hh"

namespace coot {

void add_chem_comp_sp2_C_planes_gemmi(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);

} // namespace coot