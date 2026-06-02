#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <utility>
#include <string>

namespace coot {
namespace util {

// Inline helper — same as coot::co_gemmi
inline clipper::Coord_orth co(const gemmi::Atom *at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

inline std::vector<std::pair<const gemmi::Residue *, double>>
CO_orientations_gemmi(const gemmi::Model &model) {
    std::vector<std::pair<const gemmi::Residue *, double>> scores;

    for (const auto &chain : model.chains) {
        size_t nres = chain.residues.size();
        if (nres > 2) {
            size_t res_idx_min = 1;
            size_t res_idx_max = nres - 2;
            for (size_t ires = res_idx_min; ires < res_idx_max; ++ires) {
                const gemmi::Residue &prev_res = chain.residues[ires - 1];
                const gemmi::Residue &this_res = chain.residues[ires];
                const gemmi::Residue &next_res = chain.residues[ires + 1];

                int delta_1 = this_res.seqid.num.value - prev_res.seqid.num.value;
                int delta_2 = next_res.seqid.num.value - this_res.seqid.num.value;

                if (delta_1 != 1) continue;
                if (delta_2 != 1) continue;

                const gemmi::Atom *prev_O = nullptr;
                const gemmi::Atom *prev_C = nullptr;
                const gemmi::Atom *this_O = nullptr;
                const gemmi::Atom *this_C = nullptr;
                const gemmi::Atom *next_O = nullptr;
                const gemmi::Atom *next_C = nullptr;

                for (const auto &at : prev_res.atoms) {
                    if (at.altloc == '\0') {
                        if (at.name == " C  ") prev_C = &at;
                        if (at.name == " O  ") prev_O = &at;
                    }
                }
                if (!prev_C) continue;
                if (!prev_O) continue;

                for (const auto &at : this_res.atoms) {
                    if (at.altloc == '\0') {
                        if (at.name == " C  ") this_C = &at;
                        if (at.name == " O  ") this_O = &at;
                    }
                }
                if (!this_C) continue;
                if (!this_O) continue;

                for (const auto &at : next_res.atoms) {
                    if (at.altloc == '\0') {
                        if (at.name == " C  ") next_C = &at;
                        if (at.name == " O  ") next_O = &at;
                    }
                }
                if (!next_C) continue;
                if (!next_O) continue;

                clipper::Coord_orth v1(co(prev_O) - co(prev_C));
                clipper::Coord_orth v2(co(this_O) - co(this_C));
                clipper::Coord_orth v3(co(next_O) - co(next_C));

                clipper::Coord_orth v1n(v1.unit());
                clipper::Coord_orth v2n(v2.unit());
                clipper::Coord_orth v3n(v3.unit());

                double dp_1 = clipper::Coord_orth::dot(v1n, v2n);
                double dp_2 = clipper::Coord_orth::dot(v2n, v3n);
                double sum = dp_1; // + dp_2;  (original: only dp_1)

                scores.push_back({&this_res, sum});
            }
        }
    }

    return scores;
}

} // namespace util
} // namespace coot