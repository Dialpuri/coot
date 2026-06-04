#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"

namespace coot {

struct dictionary_residue_restraints_t;

/*  gemmi port of atom_overlaps_container_t::setup_env_residue_atoms_radii

 *  Fills `neighb_atom_radius` (one value per atom in the model) with the
 *  van-der-Waals radius for each atom, looked up through the per-residue
 *  dictionary restraints.
 *
 *  If `neighb_atom_radius` is already non-empty the function returns
 *  immediately (same early-guard as the MMDB version).
 *
 *  Parameters
 *  ----------
 *  neighb_atom_radius        – out: resized to model atom count; filled
 *                              with radii (0.0 when have_dictionary == false)
 *  have_dictionary           – when false the function only resizes and
 *                              prints a message; all entries stay at 0.0
 *  neighb_dictionaries       – indexed by residue_index, used when
 *                              have_dictionary == true
 *  residue_map               – maps &atom (as unsigned long) -> residue_index
 *  type_to_vdw_radius_map    – mutable cache of type-energy → radius
 *  geom_p                    – protein geometry pointer (may be nullptr)
 *  model                     – gemmi::Model whose atoms are processed
 *
 *  Traversal order is chain → residue → atom (depth-first), matching the
 *  selection-based ordering the MMDB version produced.
 */
void atom_overlaps_container_t_setup_env_residue_atoms_radii_gemmi(
    std::vector<double>&                        neighb_atom_radius,
    bool                                        have_dictionary,
    const std::vector<dictionary_residue_restraints_t>& neighb_dictionaries,
    const std::map<unsigned long, unsigned int>& residue_map,
    std::map<std::string, double>&              type_to_vdw_radius_map,
    protein_geometry*                           geom_p,
    const gemmi::Model&                         model);

} // namespace coot