#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__is_ss_bonded/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__is_linked/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_xH_H/gemmi/function.hh"

namespace coot {
namespace reduce {

inline std::vector<gemmi::Atom *>
add_SH_H_gemmi(const std::string &H_name,
               const std::string &at_name_1,  // SG
               const std::string &at_name_2,  // CB
               const std::string &at_name_3,  // CA
               double bl,
               double angle,      // deg
               double tor_initial, // deg
               gemmi::Residue &residue,
               gemmi::CRA cra,
               const gemmi::Model &model,
               const gemmi::Structure &structure) {

    std::vector<gemmi::Atom *> empty;

    if (is_ss_bonded_gemmi(&residue, model)) {
        // don't add an H on the S
        return empty;
    }

    if (is_linked_gemmi(at_name_1, cra, structure)) {
        return empty;
    }

    return add_xH_H_gemmi(H_name, at_name_1, at_name_2, at_name_3,
                          bl, angle, tor_initial, residue);
}

} // namespace reduce
} // namespace coot