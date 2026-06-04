#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

namespace coot {

std::vector<double> get_residue_CA_position_gemmi(const std::string& cid, gemmi::Model& model) {
    std::vector<double> v;
    gemmi::Residue* residue_p = cid_to_residue_gemmi(cid, model);
    if (residue_p) {
        for (auto& atom : residue_p->atoms) {
            // gemmi does not store TER records as atoms, so no isTer() check needed
            if (atom.name == "CA") {
                v.push_back(atom.pos.x);
                v.push_back(atom.pos.y);
                v.push_back(atom.pos.z);
                break;
            }
        }
    }
    return v;
}

} // namespace coot