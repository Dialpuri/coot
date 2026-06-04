#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <iostream>
#include <string>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__reduce__position_by_bisection/gemmi/function.hh"

namespace coot {
namespace reduce {

void switch_his_protonation_gemmi(gemmi::Residue* residue_p,
                                   gemmi::Atom* current_H_atom,
                                   double bl_arom) {
    if (current_H_atom) {
        std::string atom_name = current_H_atom->name;
        std::string new_atom_name;
        if (atom_name == "HD1") new_atom_name = "HE2";
        if (atom_name == "HE2") new_atom_name = "HD1";
        if (!new_atom_name.empty()) {

            // ND1 -> HD1   and   NE2 -> HE2 by dictionary
            //
            // New atom HD1 first
            std::cout << "switch_his_protonation() " << 1 << std::endl;
            std::string at_name_1 = "CG";
            std::string at_name_2 = "ND1";
            std::string at_name_3 = "CE1";
            if (new_atom_name == "HE2") {
                at_name_1 = "CE1";
                at_name_2 = "NE2";
                at_name_3 = "CD2";
            }
            char alt_conf = current_H_atom->altloc;

            // Find atoms by name and altloc in the residue
            auto find_atom = [&](const std::string& name) -> gemmi::Atom* {
                for (auto& a : residue_p->atoms) {
                    if (a.name == name && a.altloc == alt_conf) {
                        return &a;
                    }
                }
                return nullptr;
            };

            gemmi::Atom* at_1 = find_atom(at_name_1);
            gemmi::Atom* at_2 = find_atom(at_name_2);
            gemmi::Atom* at_3 = find_atom(at_name_3);

            if (at_1 && at_2 && at_3) {
                std::cout << "switch_his_protonation() " << 2 << " " << new_atom_name << std::endl;
                current_H_atom->name = new_atom_name;
                clipper::Coord_orth pos = position_by_bisection_gemmi(at_1, at_2, at_3, bl_arom);
                double bf = current_H_atom->b_iso;
                current_H_atom->pos.x = pos.x();
                current_H_atom->pos.y = pos.y();
                current_H_atom->pos.z = pos.z();
                current_H_atom->occ  = 1.0;
                current_H_atom->b_iso = bf;
            }
        }
    }
}

} // namespace reduce
} // namespace coot