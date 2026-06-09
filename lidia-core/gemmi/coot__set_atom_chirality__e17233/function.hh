#pragma once

#include <rdkit/GraphMol/RDKitBase.h>
#include "geometry/protein-geometry.hh"
#include <string>

namespace coot {

// Set the chiral tag on an RDKit atom based on the pdbx_stereo_config
// from the dictionary atom. "R" maps to CW, "S" maps to CCW, "N" and
// absent values leave the tag untouched.
inline void set_atom_chirality_gemmi(RDKit::Atom *rdkit_at, const coot::dict_atom &dict_atom) {
   if (dict_atom.pdbx_stereo_config.first) {
      if (dict_atom.pdbx_stereo_config.second == "R") {
         RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
         rdkit_at->setChiralTag(chiral_tag);
         rdkit_at->setProp("_CIPCode", std::string("R"));
      }
      if (dict_atom.pdbx_stereo_config.second == "S") {
         RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CCW;
         rdkit_at->setProp("_CIPCode", std::string("S"));
         rdkit_at->setChiralTag(chiral_tag);
      }
      if (dict_atom.pdbx_stereo_config.second == "N") {
         // no-op: pdbx_stereo_config says N
      }
   }
}

} // namespace coot
