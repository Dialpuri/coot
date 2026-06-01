#pragma once

#include <vector>
#include <string>
#include <map>
#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "geometry/protein-geometry.hh"

class Bond_lines_container {
public:
    bool invert_deloc_bond_displacement_vector_gemmi(
        const clipper::Coord_orth &vect,
        int iat_1,
        int iat_2,
        const std::vector<const gemmi::Atom*> &residue_atoms,
        int n_atoms,
        const std::vector<coot::dict_bond_restraint_t> &bond_restraints) const;
};