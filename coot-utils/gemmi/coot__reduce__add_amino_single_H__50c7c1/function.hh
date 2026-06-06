#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_trp_indole_hydrogen/gemmi/function.hh"

namespace coot {
namespace reduce {

// 4-param overload (re-exported from existing port)
inline void add_amino_single_H_gemmi(const std::string &H_at_name,
                                     const std::string &at_name_1,
                                     const std::string &at_name_2,
                                     const std::string &at_name_3,
                                     double bl,
                                     gemmi::Residue &residue) {
    add_trp_indole_hydrogen_gemmi(H_at_name, at_name_1, at_name_2, at_name_3, bl, residue);
}

// Vector overload — the one ported here
inline void add_amino_single_H_gemmi(const std::string &H_at_name,
                                     const std::string &first_neighb,
                                     const std::vector<std::string> &second_neighb_vec,
                                     double bl,
                                     gemmi::Residue &residue) {
    if (second_neighb_vec.size() == 2) {
        add_amino_single_H_gemmi(H_at_name, second_neighb_vec[0], first_neighb, second_neighb_vec[1], bl, residue);
    }
}

} // namespace reduce
} // namespace coot