#pragma once

#include <gemmi/model.hpp>

namespace coot {

bool residues_in_order_p_gemmi(const gemmi::Chain* chain_p) {
    bool ordered_flag = true;

    if (chain_p) {
        int current_resno = -9999999;
        for (const gemmi::Residue& res : chain_p->residues) {
            int seqnum = res.seqid.num.value;
            if (seqnum < current_resno) {
                ordered_flag = false;
                break;
            } else {
                current_resno = seqnum;
            }
        }
    }
    return ordered_flag;
}

} // namespace coot