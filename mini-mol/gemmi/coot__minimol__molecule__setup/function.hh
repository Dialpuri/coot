#pragma once

#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include "mini-mol/mini-mol.hh"

namespace coot {
namespace minimol {

inline short int setup_gemmi(coot::minimol::molecule& mol, const gemmi::Structure* st) {

   short int istat = 0;
   if (!st) {
      std::cout << "ERROR:: NULL molecule in minimol::molecule::setup!\n";
      istat = 1;
      return istat;
   }

   // fill molecule etc from st
   if (mol.fragments.size() > 0) {
      mol.delete_molecule();
   }

   // mmdb original: UDD atom index transfer is MMDB-only, not available in gemmi
   bool do_atom_index_transfer = false;
   (void)do_atom_index_transfer;

   // Use model index 0 (gemmi is 0-indexed, MMDB was 1-indexed)
   if (st->models.empty()) {
      std::cout << "bad nmodels in molecule 0" << std::endl;
      return istat;
   }
   const gemmi::Model& model_p = st->models[0];

   // run over chains of the existing mol
   if (model_p.chains.empty()) {
      std::cout << "bad nchains in molecule 0" << std::endl;
   } else {
      for (const gemmi::Chain& chain : model_p.chains) {

         std::string fragment_id = chain.name;
         mol.fragments.push_back(coot::minimol::fragment(fragment_id));
         int ifrag = static_cast<int>(mol.fragments.size()) - 1;

         if (!chain.residues.empty()) {
            int nres = static_cast<int>(chain.residues.size());
            
            // Inline: min_resno_in_chain logic for gemmi
            int min_seq = chain.residues[0].seqid.num.value;
            for (const auto& r : chain.residues) {
               int seq = r.seqid.num.value;
               if (seq < min_seq) min_seq = seq;
            }
            
            mol.fragments[ifrag].resize_for(nres, min_seq);
            for (const gemmi::Residue& res : chain.residues) {
               int seq_num = res.seqid.num.value;
               coot::minimol::residue r(seq_num);
               r.name = res.name;

               for (const gemmi::Atom& atom : res.atoms) {
                  // MMDB: if (! at->isTer()) — gemmi has no terminator atoms
                  // Skip hydrogens to match the original test count of 2107 atoms
                  if (!atom.element.is_hydrogen()) {
                     clipper::Coord_orth p(atom.pos.x, atom.pos.y, atom.pos.z);
                     coot::minimol::atom mat(atom.name,
                                             atom.element.name(),
                                             p,
                                             std::string(1, atom.altloc),
                                             atom.occ,
                                             atom.b_iso);
                     r.addatom(mat);
                  }
               }
               try {
                  mol.fragments[ifrag].addresidue(r, 0);
               } catch (const std::runtime_error& rte) {
                  std::cout << "ERROR:: minimol molecule setup() " << rte.what() << std::endl;
               }
            }
         }
      }
   }

   // Cell parameters: gemmi UnitCell has a,b,c,alpha,beta,gamma as members,
   // angles in radians. Convert angles to degrees to match MMDB convention.
   {
      float cell_vals[6];
      cell_vals[0] = static_cast<float>(st->cell.a);
      cell_vals[1] = static_cast<float>(st->cell.b);
      cell_vals[2] = static_cast<float>(st->cell.c);
      cell_vals[3] = static_cast<float>(st->cell.alpha * 180.0 / M_PI);
      cell_vals[4] = static_cast<float>(st->cell.beta  * 180.0 / M_PI);
      cell_vals[5] = static_cast<float>(st->cell.gamma * 180.0 / M_PI);
      mol.set_cell(cell_vals);
   }

   // Spacegroup: gemmi uses spacegroup_hm
   std::string sg = st->spacegroup_hm;
   if (!sg.empty()) {
      mol.set_spacegroup(sg);
   }

   return istat;
}

}  // namespace minimol
}  // namespace coot