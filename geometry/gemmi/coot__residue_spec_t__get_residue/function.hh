#pragma once

#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

inline gemmi::Residue* get_residue_gemmi(const residue_spec_t& spec, gemmi::Structure* mol) {
    gemmi::Residue *r = nullptr;

    if (!mol) return r;

    if (mol->models.empty()) return r;

    auto &model = mol->models[0]; // GetModel(1) in MMDB is 1-indexed → models[0]

    for (auto &chain : model.chains) {
        std::string this_chain_id(chain.name);
        if (this_chain_id == spec.chain_id) {
            for (auto &residue : chain.residues) {
                int this_res_no = residue.seqid.num.value;
                if (this_res_no == spec.res_no) {
                    int n_atoms = static_cast<int>(residue.atoms.size());
                    if (n_atoms > 0) {
                        r = &residue;
                    }
                }
                if (r) break;
            }
        }
        if (r) break;
    }

    return r;
}

} // namespace coot