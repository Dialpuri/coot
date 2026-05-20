#pragma once

#include <vector>
#include <string>
#include "mini-mol/atom-quads.hh"

namespace coot {
namespace monomer_utils {

std::vector<atom_index_quad>
get_atom_index_quads_gemmi(const std::vector<atom_name_quad> &atom_name_quads_in,
                           const std::vector<std::string>& atom_names) {
   std::vector<atom_index_quad> v;
   for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
      int ifound = 0;
      for (size_t i1=0; i1<atom_names.size(); i1++) {
         if (atom_names[i1] == atom_name_quads_in[iquad].atom_name(0)) {
            for (size_t i2=0; i2<atom_names.size(); i2++) {
               if (atom_names[i2] == atom_name_quads_in[iquad].atom_name(1)) {
                  for (size_t i3=0; i3<atom_names.size(); i3++) {
                     if (atom_names[i3] == atom_name_quads_in[iquad].atom_name(2)) {
                        for (size_t i4=0; i4<atom_names.size(); i4++) {
                           if (atom_names[i4] == atom_name_quads_in[iquad].atom_name(3)) {
                              v.push_back(atom_index_quad(i1, i2, i3, i4));
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if (v.size() < atom_name_quads_in.size()) {
      std::cout << "Monomer utils: Failure to find correct atom quads in residue atoms\n" ;
      for (unsigned int iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
         std::cout << "  quad needed: :"
                   << atom_name_quads_in[iquad].atom_name(0) << ":  :"
                   << atom_name_quads_in[iquad].atom_name(1) << ":  :"
                   << atom_name_quads_in[iquad].atom_name(2) << ":  :"
                   << atom_name_quads_in[iquad].atom_name(3) << ":\n";
      }
      for (unsigned int iv=0; iv<v.size(); iv++) {
         std::cout << "  found quad: "
                   << v[iv].index1 << "  "
                   << v[iv].index2 << "  "
                   << v[iv].index3 << "  "
                   << v[iv].index4 << "\n";
      }
      for (size_t i1=0; i1<atom_names.size(); i1++) {
         std::cout << "  res atom " << i1 << " " << atom_names[i1] << "\n";
      }
   } else {
      // std::cout << "found all quads in residue atoms\n" ;
   } 
   return v;
}

} // namespace monomer_utils
} // namespace coot