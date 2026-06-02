#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <utility>

#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"

namespace coot {

std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>>
torsionable_bonds_monomer_internal_gemmi(
    gemmi::Residue const* residue_p,
    std::vector<gemmi::Atom*> atom_selection,
    int n_selected_atoms,
    bool include_pyranose_ring_torsions_flag,
    coot::protein_geometry *geom_p
) {
   std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> v;

   bool hydrogen_torsions = false;
   std::string rn = residue_p->name;
   std::vector<dict_torsion_restraint_t> tors_restraints =
      geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);
   bool is_pyranose = false;
   std::string group = geom_p->get_group(rn);
   // CCD dictionaries are marked at D-SACCHARIDE or SACCHARIDE (FUC)
   // SIA (sialic acid) is marked as NON-POLYMER.  Hmm..
   if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
       group == "D-SACCHARIDE" || group == "SACCHARIDE")
      is_pyranose = true;

   if (tors_restraints.size()) {
      for (unsigned int itor=0; itor<tors_restraints.size(); itor++) {

         if (!tors_restraints[itor].is_const()) {
            std::string tr_atom_name_2 = tors_restraints[itor].atom_id_2_4c();
            std::string tr_atom_name_3 = tors_restraints[itor].atom_id_3_4c();

            for (int iat1=0; iat1<n_selected_atoms; iat1++) {
               std::string atom_name_1 = coot::pad_atom_name(
                   atom_selection[iat1]->name, atom_selection[iat1]->element.name());
               for (int iat2=0; iat2<n_selected_atoms; iat2++) {
                  if (iat1 != iat2) {
                     std::string atom_name_2 = coot::pad_atom_name(
                         atom_selection[iat2]->name, atom_selection[iat2]->element.name());
                     if (atom_name_1 == tr_atom_name_2) {
                        if (atom_name_2 == tr_atom_name_3) {

                           if ((include_pyranose_ring_torsions_flag == 1) ||
                               (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
                               (!is_pyranose)) {

                              std::pair<gemmi::Atom*, gemmi::Atom*> p(atom_selection[iat1],
                                                               atom_selection[iat2]);
                              v.push_back(p);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}

} // namespace coot