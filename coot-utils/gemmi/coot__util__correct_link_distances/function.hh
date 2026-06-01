#pragma once

#include <gemmi/model.hpp>
#include <cmath>

namespace coot {
namespace util {

inline gemmi::Atom* find_atom_by_address(gemmi::Model& model, const gemmi::AtomAddress& addr) {
    for (auto& chain : model.chains) {
        if (chain.name != addr.chain_name) continue;
        for (auto& res : chain.residues) {
            if (res.seqid.num.value != addr.res_id.seqid.num.value) continue;
            for (auto& atom : res.atoms) {
                if (atom.name == addr.atom_name && atom.altloc == addr.altloc) {
                    return &atom;
                }
            }
        }
    }
    return nullptr;
}

inline void correct_link_distances_gemmi(gemmi::Structure& st) {
    for (auto& model : st.models) {
        for (auto& conn : st.connections) {
            gemmi::Atom* a1 = find_atom_by_address(model, conn.partner1);
            gemmi::Atom* a2 = find_atom_by_address(model, conn.partner2);
            if (!a1 || !a2) continue;
            gemmi::Vec3 diff = a1->pos - a2->pos;
            double real_dist = diff.length();
            if (std::fabs(conn.reported_distance - real_dist) >= 0.01) {
                conn.reported_distance = real_dist;
            }
        }
    }
}

} // namespace util
} // namespace coot