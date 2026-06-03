#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <rdkit/GraphMol/RWMol.h>

namespace coot {

// gemmi port of coot::make_residue
// Converts an RDKit molecule (with conformer) into a gemmi::Residue*.
// Returns nullptr when the molecule has no atoms.
gemmi::Residue* make_residue_gemmi(const RDKit::ROMol& rdkm, int iconf, const std::string& res_name);

} // namespace coot