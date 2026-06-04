#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <vector>
#include <set>
#include <iostream>

#include "geometry/hb-types.hh"
#include "geometry/protein-donor-acceptors.hh"

namespace coot {

void buried_unsatisfied_polar_atoms_gemmi(gemmi::Structure& st, int model_index = 0) {
   gemmi::Model& model = st.models[model_index];

   std::set<unsigned int> fixed_atom_indices; // empty

   // Collect all atoms from the model with their residue/atom names
   std::vector<const gemmi::Atom*> atom_selection;
   std::vector<std::string> res_names;
   std::vector<std::string> at_names;
   for (gemmi::Chain& chain : model.chains) {
      for (gemmi::Residue& res : chain.residues) {
         for (gemmi::Atom& at : res.atoms) {
            atom_selection.push_back(&at);
            res_names.push_back(res.name);
            at_names.push_back(at.name);
         }
      }
   }

   int n_selected_atoms = static_cast<int>(atom_selection.size());

   std::cout << "selected " << n_selected_atoms << " atoms " << std::endl;

   if (n_selected_atoms > 0) {
      // Build a contacts map using gemmi NeighborSearch (replaces contacts_by_bricks)
      // max distance between hydrogen atom and a hydrogen bond acceptor = 3.0
      gemmi::NeighborSearch ns(model, st.cell, 3.0);
      ns.populate(false);

      std::vector<std::set<unsigned int> > vec(n_selected_atoms);
      for (int i = 0; i < n_selected_atoms; i++) {
         const gemmi::Atom* at = atom_selection[i];
         std::vector<gemmi::NeighborSearch::Mark*> hits =
             ns.find_atoms(at->pos, '\0', 0.0, 3.0);
         for (gemmi::NeighborSearch::Mark* m : hits) {
            gemmi::CRA cra = m->to_cra(model);
            if (cra.atom) {
               // Find index of this atom in atom_selection
               for (int j = 0; j < n_selected_atoms; j++) {
                  if (atom_selection[j] == cra.atom) {
                     vec[i].insert(static_cast<unsigned int>(j));
                     break;
                  }
               }
            }
         }
      }

      coot::quick_protein_donor_acceptors pda;

      std::vector<bool> is_polar(n_selected_atoms, false);
      for (int i = 0; i < n_selected_atoms; i++) {
         quick_protein_donor_acceptors::key k(res_names[i], at_names[i]);
         coot::hb_t hb_type = pda.get_type(k);
         if (hb_type == HB_BOTH || hb_type == HB_DONOR || hb_type == HB_ACCEPTOR || hb_type == HB_HYDROGEN)
            is_polar[i] = true;
      }

      if (!vec.empty()) {
         int vs = static_cast<int>(vec.size());
         if (vs != n_selected_atoms) {
            std::cout << "size problem " << vs << " " << n_selected_atoms << std::endl;
         } else {
            for (int i = 0; i < n_selected_atoms; i++) {
               if (is_polar[i]) {
                  if (!vec[i].empty()) {
                     const gemmi::Atom* at = atom_selection[i];
                     if (at) {
                        bool found_something = false;
                        quick_protein_donor_acceptors::key key_1(res_names[i], at_names[i]);
                        for (unsigned int idx : vec[i]) {
                           const gemmi::Atom* at_neighb = atom_selection[idx];
                           if (at_neighb) {
                              quick_protein_donor_acceptors::key key_2(res_names[idx], at_names[idx]);
                              std::pair<bool, bool> is_valid = pda.is_hydrogen_bond_by_types(key_1, key_2);
                              if (is_valid.first)
                                 if (is_valid.second)
                                    found_something = true;
                           }
                        }
                        if (!found_something) {
                           // std::cout << "Nothing was found for atom " << res_names[i] << " " << at_names[i] << std::endl;
                        }
                     }
                  }
               }
            }
         }
      } else {
         std::cout << "empty vec - sad face" << std::endl;
      }
   }

}

} // namespace coot