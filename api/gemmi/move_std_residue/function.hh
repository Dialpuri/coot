#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <map>
#include <string>
#include <clipper/core/coords.h>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_ori_to_this_res/gemmi/function.hh"

int move_std_residue_gemmi(gemmi::CRA moving_residue, const gemmi::Residue* reference_residue);

inline int move_std_residue_gemmi(gemmi::CRA moving_residue, const gemmi::Residue* reference_residue) {
    std::map<std::string, clipper::RTop_orth> rtops =
        coot::util::get_ori_to_this_res_gemmi(*reference_residue);

    int istat = 1; // success

    if (!reference_residue) {
        istat = 0;
    } else if (rtops.empty()) {
        istat = 0;
    } else {
        if (!moving_residue.residue || moving_residue.residue->atoms.empty()) {
            istat = 0;
        } else {
            for (auto& atom : moving_residue.residue->atoms) {
                // Normalize altloc: gemmi uses ' ' or '\0' for "no altloc"
                // get_ori_to_this_res_gemmi stores keys as "" for these
                char alt_conf = atom.altloc;
                std::string alt_str;
                if (alt_conf == ' ' || alt_conf == '\0') {
                    alt_str = "";
                } else {
                    alt_str = std::string(1, alt_conf);
                }

                auto it = rtops.find(alt_str);

                if (it != rtops.end()) {
                    clipper::Coord_orth co(atom.pos.x, atom.pos.y, atom.pos.z);
                    clipper::Coord_orth rotted = co.transform(it->second);
                    atom.pos.x = rotted.x();
                    atom.pos.y = rotted.y();
                    atom.pos.z = rotted.z();
                }
            }
        }
    }

    return istat;
}