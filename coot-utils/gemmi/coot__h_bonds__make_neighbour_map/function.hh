#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

namespace coot {
namespace h_bonds {

/*
 * Build a neighbour map: for each atom in the two CRA sets, record all
 * same-residue neighbours within a covalent-bond distance range.
 *
 * First pass (cras_1):  radius 0.1..1.8  Å
 * Second pass (cras_2): radius 0.1..1.7  Å  (deduped against first pass)
 *
 * Neighbours are sorted by distance (closest first).
 */
std::map<gemmi::Atom*, std::vector<std::pair<gemmi::Atom*, float>>>
make_neighbour_map_gemmi(const std::vector<gemmi::CRA>& cras_1,
                         const std::vector<gemmi::CRA>& cras_2);

} // namespace h_bonds
} // namespace coot
