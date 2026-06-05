#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"

namespace coot {

// gemmi port of coot::restraints_container_t::add_link_torsion_for_phi_psi
//
// Scans the protein_geometry dictionary for link_type (e.g. "TRANS"),
// matches torsion restraint atom names against atoms in two consecutive
// residues (given by gemmi::CRA), and adds TORSION_RESTRAINT entries to
// restraints_vec.
//
// Returns the number of torsion restraints added.
int add_link_torsion_for_phi_psi_gemmi(
    std::vector<simple_restraint> &restraints_vec,
    const std::vector<gemmi::CRA> &cras,
    std::string link_type,
    int phi_psi_restraints_type,
    gemmi::CRA first,
    gemmi::CRA second,
    short int is_fixed_first,
    short int is_fixed_second,
    const coot::protein_geometry &geom);

}  // namespace coot