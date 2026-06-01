#pragma once

#include "geometry/protein-geometry.hh"
#include <string>
#include <vector>

namespace coot {

/**
 * Port of coot::restraints_container_t::construct_nbc_for_moving_non_moving_bonded
 *
 * Checks whether a pair of atoms (identified by name) at boundary between
 * moving and fixed regions should be added to the filtered non-bonded
 * atom indices. The pair is added unless it is already covered by a
 * bond, angle, or torsion restraint in the given link_type.
 *
 * @param filtered_non_bonded_atom_indices  output vector of vectors — jat is
 *          appended to filtered_non_bonded_atom_indices[iat] when the pair
 *          is not already constrained.
 * @param iat  index of the first atom
 * @param jat  index of the second atom
 * @param atom_name_1  name of atom iat (e.g. " N", " CA")
 * @param atom_name_2  name of atom jat
 * @param link_type    type of link (e.g. "link")
 * @param geom         protein geometry containing link definitions
 */
void construct_nbc_for_moving_non_moving_bonded_gemmi(
    std::vector<std::vector<unsigned int>>& filtered_non_bonded_atom_indices,
    unsigned int iat,
    unsigned int jat,
    const std::string& atom_name_1,
    const std::string& atom_name_2,
    const std::string& link_type,
    const coot::protein_geometry& geom);

} // namespace coot