#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <stdexcept>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_member_p/gemmi/function.hh"

namespace coot {

inline std::string chain_atoms_segid_gemmi(const gemmi::Chain* chain_p) {
    std::vector<std::string> seg_ids;

    for (const auto& residue : chain_p->residues) {
        std::string seg_id = residue.segment;
        if (seg_ids.empty()) {
            seg_ids.push_back(seg_id);
        } else {
            if (!coot::is_member_p_gemmi(seg_ids, seg_id)) {
                std::string mess = "No consistent segids for chain ";
                mess += chain_p->name;
                throw std::runtime_error(mess);
            }
        }
    }

    if (seg_ids.empty()) {
        std::string mess = "No segids for chain ";
        mess += chain_p->name;
        throw std::runtime_error(mess);
    }
    return seg_ids[0];
}

} // namespace coot