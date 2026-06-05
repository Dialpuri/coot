#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"

namespace coot {

/*
 * gemmi port of restraints_container_t::make_flanking_atoms_rama_restraints
 *
 * Original function signature:
 *   int coot::restraints_container_t::make_flanking_atoms_rama_restraints(
 *       const protein_geometry &geom)
 *
 * The original uses member variables:
 *   chain_id_save, istart_res, iend_res, istart_minus_flag, iend_plus_flag
 *
 * and calls mol->Select() to find residues, then add_rama() for each triple.
 *
 * This gemmi version takes a gemmi::Structure and all parameters explicitly.
 * Returns the number of Rama restraints added (0 in current implementation,
 * matching the original).
 */

int make_flanking_atoms_rama_restraints_gemmi(
    gemmi::Structure& mol,
    const std::string& chain_id_save,
    int istart_res,
    int iend_res,
    int istart_minus_flag,
    int iend_plus_flag,
    const coot::protein_geometry& geom);

/*
 * Helper struct that mirrors the logic: finds the triple of residues for
 * flanking Rama restraints and returns them.
 */
struct rama_flanking_triple_t {
    const gemmi::Residue* first  = nullptr;
    const gemmi::Residue* middle = nullptr;
    const gemmi::Residue* third  = nullptr;
    bool is_valid() const { return first && middle && third; }
};

std::vector<rama_flanking_triple_t>
find_rama_flanking_triples_gemmi(
    const gemmi::Model& model,
    const std::string& chain_id_save,
    int istart_res,
    int iend_res,
    int istart_minus_flag,
    int iend_plus_flag);

} // namespace coot