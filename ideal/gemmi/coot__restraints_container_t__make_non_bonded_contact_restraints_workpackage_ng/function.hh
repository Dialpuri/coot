#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <string>
#include <atomic>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__is_in_same_ring/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__check_for_O_C_1_5_relation/gemmi/function.hh"

namespace coot {

void make_non_bonded_contact_restraints_workpackage_ng_gemmi(
    int ithread,
    int imol,
    const coot::protein_geometry &geom,
    const std::vector<std::set<int> > &bonded_atom_indices,
    const std::vector<int> &angle_atom_indices,
    const std::vector<std::set<unsigned int> > &vcontacts,
    std::pair<unsigned int, unsigned int> atom_index_range_pair,
    const std::set<int> &fixed_atom_indices,
    const std::vector<std::string> &energy_type_for_atom,
    bool extended_atom_mode,
    const std::vector<gemmi::Atom *> &atom,
    const std::vector<gemmi::Residue *> &atom_residues,
    const std::vector<gemmi::Chain *> &atom_chains,
    const std::vector<int> &atom_residue_indices,
    const std::vector<bool> &atom_is_metal,
    const std::vector<bool> &atom_is_hydrogen,
    const std::vector<bool> &H_atom_parent_atom_is_donor_vec,
    const std::vector<bool> &atom_is_acceptor_vec,
    std::vector<std::set<int> > *non_bonded_contacts_atom_indices_p,
    std::vector<simple_restraint> *nbc_restraints_fragment_p,
    std::atomic<unsigned int> &done_count);

} // namespace coot
