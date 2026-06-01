#pragma once

#include <string>
#include <algorithm>
#include <gemmi/model.hpp>

namespace coot {
namespace beam_in_linked_residue {

inline void delete_atom_gemmi(gemmi::Residue *res, const std::string &atom_name) {
    if (!res) return;
    auto it = std::remove_if(
        res->atoms.begin(),
        res->atoms.end(),
        [&atom_name](const gemmi::Atom &atom) {
            return atom.name == atom_name;
        });
    if (it != res->atoms.end()) {
        res->atoms.erase(it, res->atoms.end());
    }
}

} // namespace beam_in_linked_residue
} // namespace coot