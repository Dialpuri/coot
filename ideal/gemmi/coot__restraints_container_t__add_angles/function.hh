#pragma once

#include <vector>
#include <set>
#include <string>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

int add_angles_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::set<int>& fixed_atom_indices,
    int idr,
    const std::vector<gemmi::Atom*>& res_atoms,
    gemmi::Residue* sel_res,
    const coot::protein_geometry& geom);

} // namespace coot