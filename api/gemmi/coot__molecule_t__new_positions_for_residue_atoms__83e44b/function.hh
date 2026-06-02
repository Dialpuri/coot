#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <iostream>
#include "api/moved-atom.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

namespace coot {

static std::string strip_atom_name(const std::string& name) {
    std::string trimmed = name;
    trimmed.erase(0, trimmed.find_first_not_of(' '));
    if (!trimmed.empty())
        trimmed.erase(trimmed.find_last_not_of(' ') + 1);
    return trimmed;
}

int new_positions_for_residue_atoms_gemmi(const std::string &residue_cid,
                                          const std::vector<coot::api::moved_atom_t> &moved_atoms,
                                          gemmi::Model& model) {
    gemmi::Residue* residue_p = coot::cid_to_residue_gemmi(residue_cid, model);
    if (!residue_p) {
        std::cout << "ERROR:: in new_positions_for_residue_atoms() failed to find residue " << residue_cid << std::endl;
        return -1;
    }

    if (moved_atoms.empty()) {
        return 0;
    }

    int updated = 0;
    for (const auto& moved_atom : moved_atoms) {
        std::string query_name = strip_atom_name(moved_atom.atom_name);
        for (auto& atom : residue_p->atoms) {
            if (strip_atom_name(atom.name) == query_name) {
                atom.pos.x = moved_atom.x;
                atom.pos.y = moved_atom.y;
                atom.pos.z = moved_atom.z;
                updated++;
                break;
            }
        }
    }
    return updated;
}

} // namespace coot