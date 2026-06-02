#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>
#include <map>
#include <vector>
#include <string>

#include "geometry/protein-geometry.hh"

namespace coot {

/*
 * gemmi-port of atom_overlaps_container_t::init_for_all_atom
 *
 * Populates output parameters that mirror the container's member state:
 *   - neighbours            (residue pointers)
 *   - dictionary_map        (residue-name -> restraints)
 *   - have_dictionary       (all dictionary lookups succeeded)
 *   - molecule_has_hydrogens
 *   - overlap_mode          (set to ALL_ATOM == 1)
 *   - udd_h_bond_type_handle
 *   - udd_residue_index_handle
 *
 * Uses gemmi::Model traversal instead of mmdb::Manager.
 */
void atom_overlaps_container_t_init_for_all_atom_gemmi(
    std::vector<gemmi::Residue const*>& neighbours,
    std::map<std::string, dictionary_residue_restraints_t>& dictionary_map,
    bool& have_dictionary,
    bool& molecule_has_hydrogens,
    int& overlap_mode,
    unsigned int& udd_h_bond_type_handle,
    unsigned int& udd_residue_index_handle,
    gemmi::Model const& model,
    protein_geometry* geom_p);

} // namespace coot