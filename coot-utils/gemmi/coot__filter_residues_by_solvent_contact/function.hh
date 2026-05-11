#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>

namespace coot {

// Filter residues by solvent contact - ported from MMDB to gemmi
// Returns residues that are either not HOH, or HOH within water_dist_max of any non-H/C atom in res_ref
std::vector<const gemmi::Residue*>
filter_residues_by_solvent_contact_gemmi(const gemmi::Residue& res_ref,
                                         const gemmi::Structure& st,
                                         const std::vector<const gemmi::Residue*>& residues,
                                         double water_dist_max);

} // namespace coot