#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>
#include "lidia-core/lbg-molfile.hh"
#include "geometry/protein-geometry.hh"

namespace lig_build {

// gemmi port of molfile_molecule_t constructor
// Original: lig_build::molfile_molecule_t::molfile_molecule_t(mmdb::Residue*, ...)
// Takes gemmi::Residue* instead of mmdb::Residue*
molfile_molecule_t molfile_molecule_t_gemmi(
    gemmi::Residue* residue_p,
    const coot::dictionary_residue_restraints_t &restraints);

} // namespace lig_build