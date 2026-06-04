#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>

std::pair<int, gemmi::Residue *>
find_serial_number_for_insert_gemmi(gemmi::Model &model,
                                    int seqnum_for_new,
                                    const std::string &ins_code_for_new,
                                    const std::string &chain_id) {

    int iserial_no = -1;
    gemmi::Residue *res = nullptr;

    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != chain_id)
            continue;

        int ires = 0;
        for (gemmi::Residue &residue : chain.residues) {
            int diff = residue.seqid.num.value - seqnum_for_new;

            if (diff > 0) {
                res = &residue;
                iserial_no = ires;
                break;
            } else if (diff == 0) {
                // Normalize insertion codes: MMDB "" and gemmi ' ' both mean "none"
                // Compare using normalized values
                auto norm_ins = ins_code_for_new.empty() ? " " : ins_code_for_new;
                char ins_code_this = residue.seqid.icode;
                if (ins_code_this > norm_ins[0]) {
                    res = &residue;
                    iserial_no = ires;
                    break;
                }
            }
            ++ires;
        }
    }
    return std::pair<int, gemmi::Residue *>(iserial_no, res);
}