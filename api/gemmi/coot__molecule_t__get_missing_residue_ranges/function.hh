#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include "coot/api/coot-molecule.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__closest_approach/gemmi/function.hh"

namespace coot {

// gemmi port of coot::molecule_t::get_missing_residue_ranges
// Takes a gemmi::Model and returns ranges of missing residue numbers
// between consecutive residues that are not close in space.
std::vector<residue_range_t> get_missing_residue_ranges_gemmi(const gemmi::Model& model) {
    std::vector<residue_range_t> v;

    for (const auto& chain : model.chains) {
        std::string chain_id = chain.name;
        for (size_t ires = 0; ires + 1 < chain.residues.size(); ires++) {
            const gemmi::Residue& residue_this = chain.residues[ires];
            const gemmi::Residue& residue_next = chain.residues[ires + 1];

            int rn1 = residue_this.seqid.num.value;
            int rn2 = residue_next.seqid.num.value;

            if (rn2 > (rn1 + 1)) {
                bool is_close = false;
                std::pair<bool, float> ca = closest_approach_gemmi(residue_this, residue_next);
                if (ca.first) {
                    if (ca.second < 3.0) {
                        is_close = true;
                    }
                }
                if (!is_close) {
                    int rn1_gap = rn1 + 1;
                    int rn2_gap = rn2 - 1;
                    residue_range_t rr(chain_id, rn1_gap, rn2_gap);
                    v.push_back(rr);
                }
            }
        }
    }

    return v;
}

} // namespace coot