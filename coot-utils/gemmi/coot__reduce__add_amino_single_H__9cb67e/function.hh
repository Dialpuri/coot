#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_trp_indole_hydrogen/gemmi/function.hh"

namespace coot {
namespace reduce {

inline void add_amino_single_H_gemmi(const std::string &H_at_name,
                                     const std::string &at_name_1,
                                     const std::string &at_name_2,
                                     const std::string &at_name_3,
                                     double bl,
                                     gemmi::Residue &residue) {
    add_trp_indole_hydrogen_gemmi(H_at_name, at_name_1, at_name_2, at_name_3, bl, residue);
}

} // namespace reduce
} // namespace coot