#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/contact.hpp>
#include <map>
#include <set>
#include <vector>

namespace coot {

std::map<const gemmi::Residue*, std::set<const gemmi::Residue*>>
residues_near_residues_gemmi(
    const std::vector<std::pair<bool, gemmi::CRA>>& residues_vec,
    gemmi::Structure& st,
    float dist_crit);

} // namespace coot