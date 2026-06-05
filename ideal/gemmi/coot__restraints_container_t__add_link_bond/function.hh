#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <set>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

// gemmi port of coot::restraints_container_t::add_link_bond
//
// Adds link-bond restraints between two residues based on the dictionary
// protein_geometry.  Scans all link entries in `geom` for the given
// `link_type` (e.g. "TRANS"), then for each matching link_bond_restraint
// looks for the corresponding atoms in both residues.  When a pair is
// found and their altlocs are compatible, a BOND_RESTRAINT simple_restraint
// is pushed onto `restraints_vec` and the bond is recorded in
// `bonded_atom_indices`.
//
// Returns the number of bonds added.
//
// Parameters:
//   restraints_vec         — vector to receive new simple_restraint objects
//   bonded_atom_indices    — maps each atom index to the set of atoms it is
//                            bonded to; updated in-place
//   link_type              — dictionary link type, e.g. "TRANS"
//   first                  — pointer to the first gemmi::Residue
//   second                 — pointer to the second gemmi::Residue
//   is_fixed_first         — fix flag for atoms in the first residue
//   is_fixed_second        — fix flag for atoms in the second residue
//   geom                   — protein geometry (dictionary data)
//   first_residue_index    — integer index of the first residue (used to
//                            compute global atom indices: idx = res_idx*100 +
//                            local_atom_idx)
//   second_residue_index   — same for the second residue
int add_link_bond_gemmi(
    std::vector<simple_restraint>& restraints_vec,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::string link_type,
    gemmi::Residue const *first,
    gemmi::Residue const *second,
    short int is_fixed_first,
    short int is_fixed_second,
    const coot::protein_geometry &geom,
    int first_residue_index,
    int second_residue_index);

} // namespace coot