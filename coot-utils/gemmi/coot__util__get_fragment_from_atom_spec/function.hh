#pragma once

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <optional>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__create_mmdbmanager_from_res_selection/gemmi/function.hh"

namespace coot { namespace util {

std::pair<std::optional<gemmi::Structure>, std::vector<coot::residue_spec_t>>
get_fragment_from_atom_spec_gemmi(const coot::atom_spec_t &atom_spec,
                                  const gemmi::Structure& st_in) {

   std::optional<gemmi::Structure> mol;
   std::vector<coot::residue_spec_t> v;

   // We only want one model, so take the first.
   if (st_in.models.empty()) return std::make_pair(std::nullopt, v);

   const gemmi::Model& model = st_in.models[0];

   gemmi::Residue* search_residue = nullptr;
   gemmi::Atom* search_atom = nullptr;
   gemmi::Chain* search_chain = nullptr;

   // Find the atom matching the spec
   for (const gemmi::Chain& chain : model.chains) {
      if (chain.name == atom_spec.chain_id) {
         for (const gemmi::Residue& residue : chain.residues) {
            if (residue.seqid.num.value == atom_spec.res_no) {
               for (const gemmi::Atom& atom : residue.atoms) {
                  // gemmi stores atom.name without padding (e.g. "CA" not " CA ")
                  if (atom.name == atom_spec.atom_name) {
                     // Normalize alt_conf: MMDB uses "" for no alt conf,
                     // gemmi uses '\0' for no alt conf
                     std::string alt_conf = (atom.altloc != '\0') ? std::string(1, atom.altloc) : "";
                     if (alt_conf == atom_spec.alt_conf) {
                        search_atom = const_cast<gemmi::Atom*>(&atom);
                        search_residue = const_cast<gemmi::Residue*>(&residue);
                        search_chain = const_cast<gemmi::Chain*>(&chain);
                        break;
                     }
                  }
               }
            }
            if (search_atom) break;
         }
      }
      if (search_atom) break;
   }

   if (search_atom && search_residue && search_chain) {
      int resno_this = search_residue->seqid.num.value;
      int resno_top = resno_this;
      int resno_bot = resno_this;
      gemmi::Residue* residue_top = search_residue;
      gemmi::Residue* residue_bot = search_residue;

      // Find the index of the search residue in the chain
      int search_idx = -1;
      for (int i = 0; i < static_cast<int>(search_chain->residues.size()); i++) {
         if (&search_chain->residues[i] == search_residue) {
            search_idx = i;
            break;
         }
      }

      if (search_idx >= 0) {
         // Search forwards on this chain
         for (int i = search_idx + 1; i < static_cast<int>(search_chain->residues.size()); i++) {
            int ioff = i - search_idx;
            gemmi::Residue& res = search_chain->residues[i];
            if (res.seqid.num.value == resno_this + ioff) {
               residue_top = &res;
               resno_top = res.seqid.num.value;
            } else {
               break;
            }
         }

         // Search backwards on this chain
         for (int i = search_idx - 1; i >= 0; i--) {
            int ioff = i - search_idx;
            gemmi::Residue& res = search_chain->residues[i];
            if (res.seqid.num.value == resno_this + ioff) {
               residue_bot = &res;
               resno_bot = res.seqid.num.value;
            } else {
               break;
            }
         }
      }

      if (residue_bot && residue_top) {
         // Collect all residues from bot to top (continuous range)
         std::vector<gemmi::Residue*> selected_residues;
         for (gemmi::Residue& res : search_chain->residues) {
            int rno = res.seqid.num.value;
            if (rno >= resno_bot && rno <= resno_top) {
               selected_residues.push_back(&res);
            }
         }

         if (!selected_residues.empty()) {
            auto mol_info = coot::util::create_mmdbmanager_from_res_selection_gemmi(
               st_in, selected_residues, 0, 0, "", atom_spec.chain_id, 0);

            if (mol_info.second) {
               mol = mol_info.first;
               for (gemmi::Residue* res : selected_residues) {
                  std::string icode(1, res->seqid.icode);
                  v.push_back(coot::residue_spec_t(
                     search_chain->name,
                     res->seqid.num.value,
                     icode));
               }
            }
         }

      } else {
         if (!residue_top)
            std::cout << "ERROR:: missing top residue in fragment" << std::endl;
         if (!residue_bot)
            std::cout << "ERROR:: missing bot residue in fragment" << std::endl;
      }
   }

   return std::make_pair(mol, v);
}

}} // namespace coot::util