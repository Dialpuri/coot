#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

namespace coot {

struct atom_spec_t {
    std::string atom_name;
    std::string alt_conf;
    int         res_no;
    std::string ins_code;
    std::string chain_id;
};

inline int get_asc_index_gemmi(const std::vector<gemmi::CRA>& cras,
                               const atom_spec_t& spec) {
    int asc_index = -1;
    if (spec.res_no > 0) {
        int atom_index = 0;
        for (const gemmi::CRA& cra : cras) {
            const gemmi::Atom* a = cra.atom;
            const gemmi::Residue* res = cra.residue;
            const gemmi::Chain* chain = cra.chain;
            if (res->seqid.num.value == spec.res_no &&
                (spec.ins_code.empty() ? ' ' : static_cast<char>(spec.ins_code[0])) == res->seqid.icode &&
                chain->name == spec.chain_id &&
                a->name == spec.atom_name &&
                a->altloc == static_cast<char>(spec.alt_conf.empty() ? '\0' : spec.alt_conf[0])) {
                asc_index = atom_index;
                break;
            }
            atom_index++;
        }
    }
    return asc_index;
}

}  // namespace coot