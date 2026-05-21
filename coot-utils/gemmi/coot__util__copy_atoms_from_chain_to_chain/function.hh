#pragma once
#include <gemmi/model.hpp>
#include <iostream>

namespace coot {
namespace util {

// Copy atom coordinates from one chain to another
// Both chains must have the same number of residues, and each residue
// must have the same number of atoms (in the same order)
inline void
copy_atoms_from_chain_to_chain_gemmi(gemmi::Chain *from_chain, gemmi::Chain *to_chain) {
   int n_residues_1 = from_chain->residues.size();
   int n_residues_2 = to_chain->residues.size();
   if (n_residues_2 == n_residues_1) {
      for (size_t ires = 0; ires < from_chain->residues.size(); ++ires) {
         gemmi::Residue& residue_from = from_chain->residues[ires];
         gemmi::Residue& residue_to   = to_chain->residues[ires];
         int n_atoms_1 = residue_from.atoms.size();
         int n_atoms_2 = residue_to.atoms.size();
         if (n_atoms_2 == n_atoms_1) {
            for (size_t iat = 0; iat < residue_from.atoms.size(); ++iat) {
               gemmi::Atom& at_from = residue_from.atoms[iat];
               gemmi::Atom& at_to   = residue_to.atoms[iat];
               at_to.pos = at_from.pos;
            }
         } else {
            std::cout << "ERROR:: mismatching atom count in copy_atoms_from_chain_to_chain() " << std::endl;
         }
      }
   } else {
      std::cout << "ERROR:: mismatching residue count in copy_atoms_from_chain_to_chain() " << std::endl;
   }
}

}
}