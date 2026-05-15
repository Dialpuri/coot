#pragma once

#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
namespace lidia_utils {

inline gemmi::Residue *get_residue_gemmi(const residue_spec_t &res_spec,
                                         gemmi::Structure *st) {
    if (!st) return nullptr;

    if (st->models.empty()) return nullptr;
    auto &model = st->models[0];  // GetModel(1) in MMDB is 1-indexed → models[0]

    for (auto &chain : model.chains) {
        std::string mol_chain(chain.name);
        if (mol_chain == res_spec.chain_id) {
            for (auto &residue : chain.residues) {
                if (residue.seqid.num.value == res_spec.res_no) {
                    // MMDB uses "" for no ins code; gemmi uses ' '
                    std::string gemmi_ic(1, residue.seqid.icode);
                    std::string norm_query = res_spec.ins_code.empty()
                                                 ? std::string(" ")
                                                 : res_spec.ins_code;
                    if (norm_query == gemmi_ic) {
                        return &residue;
                    }
                }
            }
        }
    }
    return nullptr;
}

} // namespace lidia_utils
} // namespace coot