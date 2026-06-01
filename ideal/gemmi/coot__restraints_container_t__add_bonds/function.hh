#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "geometry/protein-geometry.hh"
#include "geometry/hb-types.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

// add_bonds_gemmi — gemmi port of coot::restraints_container_t::add_bonds()
//
// Scans the dictionary bond restraints for residue `sel_res` (given by
// dictionary index `idr` in `geom`) against the atoms in `res_atoms`.
// For each matched bond pair that passes the altloc check, a
// simple_restraint is pushed into `restraints_vec`.
//
// Returns the number of bond restraints added.
int add_bonds_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::map<gemmi::Atom*, hb_t>& H_atom_parent_energy_type_atom_map,
    int idr,
    const std::vector<gemmi::Atom*>& res_atoms,
    gemmi::Residue* sel_res,
    const coot::protein_geometry& geom);

} // namespace coot