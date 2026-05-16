#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace typed_distances {

enum atom_type_t { NONE, C, O, N };

inline atom_type_t get_type_gemmi(gemmi::Atom* at) {
    atom_type_t t1(NONE);
    if (at == nullptr) return t1;
    std::string ele = at->element.name();
    if (ele == "C") t1 = atom_type_t(C);
    if (ele == "O") t1 = atom_type_t(O);
    if (ele == "S") t1 = atom_type_t(O);
    if (ele == "N") t1 = atom_type_t(N);
    return t1;
}

}  // namespace typed_distances
}  // namespace coot