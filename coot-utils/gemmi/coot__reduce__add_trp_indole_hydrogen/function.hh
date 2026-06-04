#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include <clipper/core/coords.h>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_alt_confs/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bisection/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__add_hydrogen_atom/gemmi/function.hh"

namespace coot {
namespace reduce {

inline const gemmi::Atom* find_atom_in_residue(const gemmi::Residue& res,
                                                const std::string& atom_name,
                                                const std::string& altconf) {
    std::string trimmed_name = trim(atom_name);
    for (const auto& atom : res.atoms) {
        if (trim(atom.name) != trimmed_name) continue;
        if (!altconf.empty()) {
            if (atom.altloc == altconf[0]) return &atom;
        } else {
            if (atom.altloc == ' ' || atom.altloc == '\0') return &atom;
        }
    }
    return nullptr;
}

void add_trp_indole_hydrogen_gemmi(const std::string& H_name,
                                    const std::string& at_name_1,
                                    const std::string& at_name_2,
                                    const std::string& at_name_3,
                                    double bl,
                                    gemmi::Residue& residue) {
    std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs_gemmi(residue);
    for (size_t i = 0; i < alt_confs.size(); i++) {
        const gemmi::Atom* at_1 = find_atom_in_residue(residue, at_name_1, alt_confs[i]);
        const gemmi::Atom* at_2 = find_atom_in_residue(residue, at_name_2, alt_confs[i]);
        const gemmi::Atom* at_3 = find_atom_in_residue(residue, at_name_3, alt_confs[i]);
        if (at_1 && at_2 && at_3) {
            clipper::Coord_orth H_pos = coot::reduce::position_by_bisection_gemmi(at_1, at_2, at_3, bl);
            float bf = static_cast<float>(at_2->b_iso);
            coot::reduce::add_hydrogen_atom_gemmi(H_name, H_pos, bf, alt_confs[i], residue);
        }
    }
}

} // namespace reduce
} // namespace coot