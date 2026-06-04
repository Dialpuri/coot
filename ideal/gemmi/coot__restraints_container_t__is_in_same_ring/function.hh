#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

bool is_in_same_ring_gemmi(
    int imol,
    gemmi::Residue const& residue,
    std::map<std::string, std::pair<bool, std::vector<std::list<std::string>>>>& residue_ring_map_cache,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const coot::protein_geometry& geom);

} // namespace coot