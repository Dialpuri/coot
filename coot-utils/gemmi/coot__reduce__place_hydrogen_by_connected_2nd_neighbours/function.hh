#pragma once

#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include <vector>

namespace coot::reduce {

inline void place_hydrogen_by_connected_2nd_neighbours_gemmi(
    unsigned int iat,
    unsigned int iat_neighb,
    const coot::dictionary_residue_restraints_t &rest,
    gemmi::Residue &residue) {

    (void)iat;
    (void)residue;

    std::vector<unsigned int> neighbs = rest.neighbours(iat_neighb, false);
    std::string ele_neighb = rest.atom_info[iat_neighb].type_symbol;
    (void)ele_neighb;

    // something complicated
    for (unsigned int i = 0; i < neighbs.size(); i++) {
    }
}

}  // namespace coot::reduce
