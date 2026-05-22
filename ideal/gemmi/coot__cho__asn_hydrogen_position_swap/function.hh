#pragma once

#include <vector>
#include <utility>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot { namespace cho {

// Port of coot::cho::asn_hydrogen_position_swap from MMDB to gemmi
// Swaps HD21 and HD22 atom positions between ASN and NAG residues
inline void asn_hydrogen_position_swap_gemmi(std::vector<std::pair<bool, gemmi::Residue*>> residues) {
    if (residues[0].second) {
        if (residues[1].second) {
            std::string rn0 = residues[0].second->name;
            std::string rn1 = residues[1].second->name;
            gemmi::Residue* r_0 = nullptr;
            gemmi::Residue* r_1 = nullptr;
            if (rn0 == "ASN") {
                if (rn1 == "NAG") {
                    r_0 = residues[0].second;
                    r_1 = residues[1].second;
                }
            }
            if (rn1 == "ASN") {
                if (rn0 == "NAG") {
                    r_1 = residues[0].second;
                    r_0 = residues[1].second;
                }
            }

            if (r_1 && r_0) {
                gemmi::Atom* at_hd21 = nullptr;
                gemmi::Atom* at_hd22 = nullptr;
                for (gemmi::Atom& atom : r_0->atoms) {
                    if (atom.name == "HD21") at_hd21 = &atom;
                    if (atom.name == "HD22") at_hd22 = &atom;
                }
                if (at_hd21 && at_hd22) {
                    clipper::Coord_orth co21 = coot::co_gemmi(at_hd21);
                    clipper::Coord_orth co22 = coot::co_gemmi(at_hd22);
                    at_hd21->pos.x = co22.x();
                    at_hd21->pos.y = co22.y();
                    at_hd21->pos.z = co22.z();
                    at_hd22->pos.x = co21.x(); // this atom will be deleted.
                    at_hd22->pos.y = co21.y();
                    at_hd22->pos.z = co21.z();
                }
            }
        }
    }
}

}} // namespace coot::cho