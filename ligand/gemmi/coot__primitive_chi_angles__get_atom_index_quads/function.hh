#pragma once
#include <string>
#include <vector>

namespace coot {

struct atom_name_quad {
   std::string names[4];
   atom_name_quad(const std::string& n0, const std::string& n1, const std::string& n2, const std::string& n3) {
      names[0] = n0; names[1] = n1; names[2] = n2; names[3] = n3;
   }
   std::string atom_name(int i) const { return names[i]; }
};

struct atom_index_quad {
   int index1, index2, index3, index4;
   atom_index_quad(int i1, int i2, int i3, int i4)
      : index1(i1), index2(i2), index3(i3), index4(i4) {}
};

namespace primitive_chi_angles {

std::vector<atom_index_quad>
get_atom_index_quads_gemmi(const std::vector<atom_name_quad> &atom_name_quads_in,
                          const std::vector<std::string> &atom_names) {
   std::vector<atom_index_quad> v;
   for (size_t iquad=0; iquad<atom_name_quads_in.size(); iquad++) {
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
   return v;
}

} // namespace primitive_chi_angles
} // namespace coot