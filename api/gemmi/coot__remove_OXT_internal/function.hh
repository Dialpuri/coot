#pragma once
#include <gemmi/model.hpp>
#include <string>

namespace coot {

inline void remove_OXT_internal_gemmi(gemmi::Residue* residue) {
    for (auto it = residue->atoms.begin(); it != residue->atoms.end(); ++it) {
        std::string name = it->name;
        // Strip trailing spaces to match standard atom name format
        while (!name.empty() && name.back() == ' ')
            name.pop_back();
        if (name == "OXT") {
            residue->atoms.erase(it);
            break;
        }
    }
}

}