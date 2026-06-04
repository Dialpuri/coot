#pragma once

#include <gemmi/model.hpp>
#include <coot/geometry/residue-and-atom-specs.hh>
#include <coot/ideal/simple-restraint.hh>
#include <vector>
#include <set>

// Carry parent context alongside each gemmi::Atom* (gemmi has no parent pointers)
struct atom_entry {
    gemmi::Atom* atom;
    const gemmi::Residue* residue;
    const gemmi::Chain* chain;
};

namespace coot {

// Build an atom_spec_t from gemmi Atom + Residue + Chain (gemmi has no parent pointers)
coot::atom_spec_t atom_spec_from_gemmi(const gemmi::Atom& atom,
                                       const gemmi::Residue& res,
                                       const gemmi::Chain& chain,
                                       int model_number = 1);

// Port of coot::restraints_container_t::add_header_metal_link_bond_ng
// Returns true if the bond was added, false if either atom was not found.
bool add_header_metal_link_bond_ng_gemmi(
    const std::vector<atom_entry>& atoms,
    std::vector<std::set<int>>& bonded_atom_indices,
    std::vector<simple_restraint>& restraints_vec,
    const coot::atom_spec_t& atom_spec_1,
    const coot::atom_spec_t& atom_spec_2,
    double dist);

} // namespace coot