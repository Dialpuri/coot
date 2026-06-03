#pragma once

#include <gemmi/model.hpp>
#include <gemmi/unitcell.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include "clipper/core/coords.h"

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_cell_symm/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__symmetry_move_atoms/gemmi/function.hh"

namespace coot {
namespace util {

int move_waters_around_protein_gemmi(gemmi::Structure& st) {

   int n_moved = 0;
   std::vector<clipper::Coord_orth> protein_coords;
   std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>> water_atoms;

   // Collect protein and water atoms from model 0
   if (!st.models.empty()) {
      gemmi::Model& model = st.models[0];
      for (gemmi::Chain& chain : model.chains) {
         for (gemmi::Residue& residue : chain.residues) {
            std::string residue_name(residue.name);
            if (residue_name == "WAT" || residue_name == "HOH") {
               for (gemmi::Atom& atom : residue.atoms) {
                  // In gemmi, TER records are not parsed as atoms, so no isTer check needed
                  clipper::Coord_orth c(atom.pos.x, atom.pos.y, atom.pos.z);
                  water_atoms.push_back(std::make_pair(&atom, c));
               }
            } else {
               for (gemmi::Atom& atom : residue.atoms) {
                  std::string ele(atom.element.name());
                  if (ele != "C") {
                     clipper::Coord_orth pt(atom.pos.x, atom.pos.y, atom.pos.z);
                     protein_coords.push_back(pt);
                  }
               }
            }
         }
      }
   }

   // OK, so water_atoms and protein atoms are filled.

   try {
      std::pair<clipper::Cell, clipper::Spacegroup> csp = get_cell_symm_gemmi(st);
      clipper::Cell cell = csp.first;
      clipper::Spacegroup spacegroup = csp.second;

      if (!cell.is_null() && !spacegroup.is_null()) {
         std::vector<std::pair<gemmi::Atom*, clipper::Coord_orth>> water_atoms_moved =
            symmetry_move_atoms_gemmi(protein_coords, water_atoms, cell, spacegroup);

         for (const auto& wp : water_atoms_moved) {
            if (wp.first) {
               wp.first->pos.x = wp.second.x();
               wp.first->pos.y = wp.second.y();
               wp.first->pos.z = wp.second.z();
               n_moved++;
            }
         }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   }

   return n_moved;
}

} // namespace util
} // namespace coot