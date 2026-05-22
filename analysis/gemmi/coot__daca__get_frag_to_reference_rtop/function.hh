#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <string>
#include <vector>
#include <map>

namespace coot { namespace daca {

// Port of get_frag_to_reference_rtop from MMDB to gemmi
// Takes fragment atoms in gemmi format and returns a rigid-body transformation
// to the reference fragment, or (false, RTop_orth()) on error.
inline std::pair<bool, clipper::RTop_orth>
get_frag_to_reference_rtop_gemmi(const std::string &res_name,
                                 const unsigned int &frag_idx,
                                 const std::vector<const gemmi::Atom*> &fragment_atoms,
                                 const std::map<std::string, std::vector<std::vector<clipper::Coord_orth>>> &reference_fragments) {
    clipper::RTop_orth rtop;
    bool status = true;
    auto it = reference_fragments.find(res_name);
    if (it != reference_fragments.end()) {
        if (frag_idx < it->second.size()) {
            const std::vector<clipper::Coord_orth> &ref_atom_positions = it->second[frag_idx];
            // convert fragment_atoms to coordinates
            std::vector<clipper::Coord_orth> residue_fragment_atoms;
            for (unsigned int i = 0; i < fragment_atoms.size(); i++) {
                const gemmi::Atom* at = fragment_atoms[i];
                clipper::Coord_orth pos(at->pos.x, at->pos.y, at->pos.z);
                residue_fragment_atoms.push_back(pos);
            }
            if (ref_atom_positions.size() == residue_fragment_atoms.size()) {
                clipper::RTop_orth rtop_1(residue_fragment_atoms, ref_atom_positions);
                rtop = rtop_1;
            } else {
                status = false;
            }
        } else {
            status = false;
        }
    } else {
        status = false;
    }
    return std::make_pair(status, rtop);
}

}} // namespace coot::daca