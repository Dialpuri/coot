#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

namespace coot {

/**
 * Port of restraints_container_t::make_strand_pseudo_bond_restraints to gemmi
 * 
 * This method creates pseudo-bond restraints for strand structures:
 * - O-O bonds (distance ~4.64 Å) 
 * - O-O-O angles (angle ~98°)
 * - CA-CA-CA angles (angle ~120°)
 * 
 * The algorithm requires consecutive residues (seqNum differs by 1).
 * 
 * @param residues Vector of residues to process (must be consecutive in seqNum)
 */
void make_strand_pseudo_bond_restraints_gemmi(
    const std::vector<const gemmi::Residue*>& residues);

} // namespace coot