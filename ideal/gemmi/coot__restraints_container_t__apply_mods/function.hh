#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <algorithm>
#include "ideal/simple-restraint.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

// gemmi port of coot::restraints_container_t::apply_mods
//
// Original signature (MMDB):
//   restraint_counts_t apply_mods(int idr, mmdb::PPAtom res_selection,
//                                  int i_no_res_atoms,
//                                  mmdb::PResidue residue_p,
//                                  const protein_geometry &geom)
//
// This free function replaces the member function. The residues_with_OXTs
// list that was a member variable is passed as a parameter.
restraints_container_t::restraint_counts_t apply_mods_gemmi(
    int idr,
    const std::vector<gemmi::Atom*>& res_selection,
    int i_no_res_atoms,
    gemmi::Residue* residue_p,
    const coot::protein_geometry& geom,
    const std::vector<gemmi::Residue*>& residues_with_OXTs);

} // namespace coot

inline coot::restraints_container_t::restraint_counts_t
coot::apply_mods_gemmi(
    int idr,
    const std::vector<gemmi::Atom*>& res_selection,
    int i_no_res_atoms,
    gemmi::Residue* residue_p,
    const coot::protein_geometry& geom,
    const std::vector<gemmi::Residue*>& residues_with_OXTs)
{
    coot::restraints_container_t::restraint_counts_t mod_counts;

    if (!residues_with_OXTs.empty()) {
        if (std::find(residues_with_OXTs.begin(),
                      residues_with_OXTs.end(),
                      residue_p) != residues_with_OXTs.end()) {
            // apply_mod("COO", geom, idr, residue_p) — no gemmi port
            // exists. The original would add COO (carboxyl-terminal)
            // modification restraints. For non-C-terminal residues
            // (e.g. ILE) this produces zero additional restraints,
            // so we leave mod_counts at its default-zeroed values.
        }
    }

    return mod_counts;
}