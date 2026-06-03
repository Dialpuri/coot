#pragma once

#include <gemmi/model.hpp>
#include <string>

namespace coot {
namespace crankshaft {

inline gemmi::Atom *get_atom_gemmi(gemmi::Residue *res_1, const std::string &atom_name_in) {
    gemmi::Atom *r = nullptr;
    // Trim whitespace from atom_name_in to match gemmi's unpadded atom.name
    auto trim = [](const std::string &s) {
        size_t start = s.find_first_not_of(' ');
        if (start == std::string::npos) return std::string();
        size_t end = s.find_last_not_of(' ');
        return s.substr(start, end - start + 1);
    };
    std::string trimmed = trim(atom_name_in);
    for (auto &at : res_1->atoms) {
        if (at.name == trimmed) {
            r = &at;
            break;
        }
    }
    return r;
}

} // namespace crankshaft
} // namespace coot