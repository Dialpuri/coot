#pragma once

#include <vector>
#include <utility>
#include <stdexcept>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__dipole__init/gemmi/function.hh"

namespace coot {
namespace dipole {

inline clipper::Coord_orth dipole_gemmi(
    std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::CRA>> dict_res_pairs
) {
    dipole_state_t state;
    init_gemmi(state, dict_res_pairs);
    if (!state.dipole_is_good_flag) {
        throw std::runtime_error("dipole not computed");
    }
    return state.dipole_;
}

} // namespace dipole
} // namespace coot
