#pragma once

#include <gemmi/model.hpp>
#include <algorithm>
#include <iostream>

namespace coot {

void sort_residues_gemmi(gemmi::Structure *mol) {

   if (mol) { 
      for (gemmi::Model &model_p : mol->models) {
         // run over chains of the existing mol
         for (gemmi::Chain &chain_p : model_p.chains) {
            // Sort residues by sequence number (ascending), then by insertion code (descending)
            std::sort(chain_p.residues.begin(), chain_p.residues.end(),
                [](const gemmi::Residue &a, const gemmi::Residue &b) {
                    int r1 = a.seqid.num.value;
                    int r2 = b.seqid.num.value;
                    if (r1 != r2) return r1 < r2;
                    // Same seqnum: compare insertion code (descending order)
                    return a.seqid.icode > b.seqid.icode;
                });
         }
      }
   } else {
      std::cout << "ERROR:: (trapped) sort_residues() called with Null molecule"
                << std::endl;
   } 
}

} // namespace coot