#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include "geometry/protein-geometry.hh"

struct gemmi_plane_restraint_data_t {
   std::vector<std::pair<int, double> > pos_sigma;
   std::vector<bool> fixed_flags;
};

namespace coot {

int add_link_plane_gemmi(
    std::string link_type,
    gemmi::CRA first,
    gemmi::CRA second,
    short int is_fixed_first_res,
    short int is_fixed_second_res,
    const coot::protein_geometry& geom,
    std::vector<gemmi_plane_restraint_data_t>& planes_out);

} // namespace coot