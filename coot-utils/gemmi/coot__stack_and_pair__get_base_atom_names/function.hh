#pragma once

#include <string>
#include <vector>
#include <set>
#include <gemmi/model.hpp>

namespace coot {
namespace stack_and_pair {

// Base atom names for nucleobases (from the original MMDB function).
static const std::set<std::string> base_atom_name_set = {
    " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N7 ", " C8 ", " N9 "
};

std::vector<std::string>
get_base_atom_names_gemmi(const gemmi::Residue &residue) {
    std::vector<std::string> v;
    v.reserve(6);
    for (const gemmi::Atom &at : residue.atoms) {
        if (base_atom_name_set.find(at.name) != base_atom_name_set.end()) {
            v.push_back(at.name);
        }
    }
    return v;
}

} // namespace stack_and_pair
} // namespace coot