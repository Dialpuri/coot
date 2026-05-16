#pragma once

#include <vector>
#include <string>
#include "mini-mol/atom-quads.hh"

namespace coot {
namespace monomer_utils {

inline std::vector<atom_index_quad>
get_atom_index_quads_gemmi(const std::vector<atom_name_quad>& atom_name_quads_in,
                           const std::vector<std::string>& atom_names) {
    std::vector<atom_index_quad> v;
    for (const auto& quad : atom_name_quads_in) {
        for (int i1 = 0; i1 < (int)atom_names.size(); i1++) {
            if (atom_names[i1] == quad.atom_name(0)) {
                for (int i2 = 0; i2 < (int)atom_names.size(); i2++) {
                    if (atom_names[i2] == quad.atom_name(1)) {
                        for (int i3 = 0; i3 < (int)atom_names.size(); i3++) {
                            if (atom_names[i3] == quad.atom_name(2)) {
                                for (int i4 = 0; i4 < (int)atom_names.size(); i4++) {
                                    if (atom_names[i4] == quad.atom_name(3)) {
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

} // namespace monomer_utils
} // namespace coot