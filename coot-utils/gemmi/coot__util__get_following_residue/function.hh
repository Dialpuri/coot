#pragma once

#include <gemmi/model.hpp>
#include <optional>
#include "geometry/residue-and-atom-specs.hh"

namespace coot::util {

// Returns a gemmi::CRA (Chain*, Residue*, Atom*) for the residue
// following the one specified by `rs` in `st`.
// atom is nullptr since the original only needs the Residue pointer.
inline std::optional<gemmi::CRA> get_following_residue_gemmi(
        const residue_spec_t &rs, gemmi::Structure &st) {
    if (st.models.empty()) return std::nullopt;

    gemmi::Model &model = st.models[0]; // model 1 in MMDB → models[0] gemmi
    bool found_this_res = false;

    for (gemmi::Chain &chain : model.chains) {
        if (chain.name != rs.chain_id) continue;
        for (size_t ires = 0; ires < chain.residues.size(); ++ires) {
            gemmi::Residue &residue = chain.residues[ires];
            if (!found_this_res) {
                if (rs.res_no == residue.seqid.num.value) {
                    // Check insertion code if present in the spec
                    if (rs.ins_code.empty() ||
                        rs.ins_code == std::string(1, residue.seqid.icode)) {
                        found_this_res = true;
                    }
                }
            } else {
                // We already found the target — this is the following residue
                return gemmi::CRA{&chain, &residue, nullptr};
            }
        }
    }
    return std::nullopt;
}

} // namespace coot::util