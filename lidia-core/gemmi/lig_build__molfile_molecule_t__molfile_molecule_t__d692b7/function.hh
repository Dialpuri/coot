#pragma once

#include "lidia-core/lbg-molfile.hh"
#include "geometry/protein-geometry.hh"
#include "clipper/core/coords.h"

#include <map>

namespace lig_build {

// gemmi-ported constructor: builds molfile_molecule_t from dictionary restraints
inline molfile_molecule_t molfile_molecule_t_gemmi(const coot::dictionary_residue_restraints_t &restraints) {

   molfile_molecule_t mol;

   std::map<std::string, int> atom_name_index;

   molfile_atom_t blank_atom(0, 0, 0, "");
   mol.atoms.push_back(blank_atom); // blank atom for 0-index.

   for (unsigned int iat = 0; iat < restraints.atom_info.size(); iat++) {
      molfile_atom_t atom(clipper::Coord_orth(0, 0, 0),
                          restraints.atom_info[iat].type_symbol,
                          restraints.atom_info[iat].atom_id);
      mol.atoms.push_back(atom);
      atom_name_index[restraints.atom_info[iat].atom_id] = static_cast<int>(iat);
   }

   for (unsigned int ibond = 0; ibond < restraints.bond_restraint.size(); ibond++) {
      const coot::dict_bond_restraint_t &bond_restraint = restraints.bond_restraint[ibond];
      auto it_1_atom_name_index = atom_name_index.find(bond_restraint.atom_id_1_4c());
      auto it_2_atom_name_index = atom_name_index.find(bond_restraint.atom_id_2_4c());
      if (it_1_atom_name_index != atom_name_index.end()) {
         if (it_2_atom_name_index != atom_name_index.end()) {
            int idx_1 = it_1_atom_name_index->second;
            int idx_2 = it_2_atom_name_index->second;

            // Inline get_bond_type (private member of molfile_molecule_t)
            bond_t::bond_type_t bond_type = lig_build::bond_t::BOND_UNDEFINED;
            std::string btype = bond_restraint.type();
            if (btype == "single")
               bond_type = lig_build::bond_t::SINGLE_BOND;
            else if (btype == "double")
               bond_type = lig_build::bond_t::DOUBLE_BOND;
            else if (btype == "triple")
               bond_type = lig_build::bond_t::TRIPLE_BOND;
            else if (btype == "aromatic")
               bond_type = lig_build::bond_t::AROMATIC_BOND;

            molfile_bond_t bond(idx_1 + 1, idx_2 + 1, bond_type);
            mol.bonds.push_back(bond);
         }
      }
   }

   return mol;
}

} // namespace lig_build