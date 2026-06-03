#pragma once

#include <map>
#include <set>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-dev/coot/geometry/residue-and-atom-specs-gemmi.hh"

namespace coot {

std::map<coot::residue_spec_t, std::set<coot::residue_spec_t>>
residues_near_residues_gemmi(gemmi::Structure* st, float dist_crit);

} // namespace coot