#pragma once

#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__dipole__charged_atoms__7dad65/gemmi/function.hh"

namespace coot {
namespace dipole {

struct dipole_state_t {
    bool dipole_is_good_flag;
    clipper::Coord_orth dipole_;
    clipper::Coord_orth residue_centre;
};

void init_gemmi(
    dipole_state_t& result,
    std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> dict_res_pairs
);

} // namespace dipole
} // namespace coot