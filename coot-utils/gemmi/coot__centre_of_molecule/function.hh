#pragma once
#include <utility>
#include <clipper/core/coords.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>

namespace coot {

std::pair<bool, clipper::Coord_orth>
centre_of_molecule_gemmi(const gemmi::Structure& st) {
   bool status = false;
   clipper::Coord_orth centre(0,0,0);
   double xs=0, ys=0, zs=0;
   int n_atoms = 0;

   for (const gemmi::Model& model : st.models) {
      for (const gemmi::Chain& chain : model.chains) {
         for (const gemmi::Residue& residue : chain.residues) {
            for (const gemmi::Atom& atom : residue.atoms) {
               // Skip hydrogens like the original skips isTer() atoms
               if (!atom.element.is_hydrogen()) {
                  xs += atom.pos.x;
                  ys += atom.pos.y;
                  zs += atom.pos.z;
                  n_atoms++;
               }
            }
         }
      }
   }

   if (n_atoms > 0) {
      status = true;
      double dna = static_cast<double>(n_atoms);
      centre = clipper::Coord_orth(xs/dna, ys/dna, zs/dna);
   }

   return std::make_pair(status, centre);
}

} // namespace coot