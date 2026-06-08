#pragma once

#include <gemmi/cifdoc.hpp>
#include <coot/geometry/protein-geometry.hh>
#include <map>
#include <vector>

namespace coot {

void add_chem_links_gemmi(
    gemmi::cif::Loop& loop,
    std::map<unsigned int, std::vector<chem_link>>& chem_link_map);

} // namespace coot
