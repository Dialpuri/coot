#pragma once

#include <gemmi/model.hpp>
#include <algorithm>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"

namespace coot {
namespace util {

// Copy a residue and delete all hydrogen atoms from the copy.
// Returns a newly allocated gemmi::Residue (caller owns with delete).
gemmi::Residue* copy_and_delete_hydrogens_gemmi(const gemmi::Residue* residue_in) {
    gemmi::Residue* copy = deep_copy_this_residue_gemmi(residue_in);
    if (!copy) {
        return nullptr;
    }

    // Remove atoms whose element is H or D (deuterium)
    copy->atoms.erase(
        std::remove_if(copy->atoms.begin(), copy->atoms.end(),
            [](const gemmi::Atom& a) {
                std::string nm = a.element.name();
                return nm == "H" || nm == "D";
            }),
        copy->atoms.end()
    );

    return copy;
}

} // namespace util
} // namespace coot