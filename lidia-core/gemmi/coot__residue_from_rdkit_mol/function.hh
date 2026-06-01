#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/unitcell.hpp>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <string>
#include <vector>
#include <memory>

namespace coot {

std::unique_ptr<gemmi::Residue> residue_from_rdkit_mol_gemmi(
    const RDKit::ROMol &mol, int conf_id, const std::string &new_comp_id);

} // namespace coot