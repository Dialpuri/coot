#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <cmath>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__get_chiral_hydrogen_index/gemmi/function.hh"

namespace coot {

int add_chirals_gemmi(
    int idr,
    std::vector<gemmi::Atom *> res_selection,
    int i_no_res_atoms,
    gemmi::CRA cra,
    const coot::protein_geometry &geom,
    std::vector<coot::simple_restraint> &restraints_vec
);

} // namespace coot