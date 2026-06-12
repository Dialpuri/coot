#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "geometry/protein-geometry.hh"
#include <iosfwd>

namespace pli {

class pi_stacking_instance_t {
public:

   // CATION_PI_STACKING sets ligand_cationic_atom_name, not the
   // ligand_ring_atom_names vector.
   //
   enum stacking_t {
     NO_STACKING,
     PI_PI_STACKING,
     PI_CATION_STACKING, // for cations on the protein residues (ligand pi)
     CATION_PI_STACKING, // for cations on the ligand (protein TRY, PRO, TRP)
   };
   gemmi::Residue *res;
   stacking_t type; // pi-pi or pi-cation
   std::vector<std::string> ligand_ring_atom_names;
   float overlap_score;
   std::string ligand_cationic_atom_name; // for cations on the ligand

   pi_stacking_instance_t(gemmi::Residue *res_in, stacking_t type_in,
                   const std::vector<std::string> &ring_atoms) : ligand_ring_atom_names(ring_atoms) {
     res = res_in;
     type = type_in;
     overlap_score = 0;
   }

   // and the constructor for CATION_PI_STACKING
   //
   pi_stacking_instance_t(gemmi::Residue *residue_in,
                   const std::string &ligand_atom_name_in) : ligand_cationic_atom_name(ligand_atom_name_in) {
     type = CATION_PI_STACKING;
     res = residue_in;
     overlap_score = 0;
   }
   friend std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);
};
std::ostream& operator<< (std::ostream& s, const pi_stacking_instance_t &spec);

class pi_stacking_container_t {
public:
   std::vector<pi_stacking_instance_t> stackings;

   size_t size() const { return stackings.size(); }
};

/*
 * gemmi-port of pli::pi_stacking_container_t constructor.
 *
 * Original:
 *   pi_stacking_container_t(const coot::dictionary_residue_restraints_t &,
 *                           const std::vector<mmdb::Residue*>,
 *                           mmdb::Residue*)
 *
 * gemmi version uses gemmi::Residue* for residue pointers.
 */
pli::pi_stacking_container_t pi_stacking_container_t_gemmi(
    const coot::dictionary_residue_restraints_t &monomer_restraints,
    const std::vector<gemmi::Residue*> &filtered_residues,
    gemmi::Residue *res_ref,
    const std::vector<std::vector<std::string>> &aromatic_ring_list);

} // namespace pli
