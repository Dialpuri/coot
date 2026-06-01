#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

class b_factor_histogram {
public:
    int n_bins;
    int n_atoms;
    float b_max;
    std::vector<std::vector<float>> b_vector;

    // chain_id == 0 means all chains
    static b_factor_histogram b_factor_histogram_gemmi(
        const gemmi::Model& model,
        char chain_id);
};

inline b_factor_histogram b_factor_histogram::b_factor_histogram_gemmi(
    const gemmi::Model& model,
    char chain_id)
{
    b_factor_histogram hist;
    // init() equivalent
    hist.n_bins  = 0;
    hist.n_atoms = 0;
    hist.b_max   = -1.0f;

    // ---- first pass: count atoms with b >= 0 and find b_max ----
    for (const auto& chain : model.chains) {
        if (chain_id != 0 && chain.name[0] != chain_id)
            continue;
        for (const auto& res : chain.residues) {
            for (const auto& atom : res.atoms) {
                float b = atom.b_iso;
                if (b >= 0.0f) {
                    hist.n_atoms++;
                    if (b > hist.b_max)
                        hist.b_max = b;
                }
            }
        }
    }

    // ---- compute n_bins ----
    // get_n_bins() always returns 40
    if (hist.n_atoms > 0) {
        hist.n_bins = 40;
    }
    hist.b_vector.resize(hist.n_bins);

    // ---- second pass: distribute b-factors into bins ----
    for (const auto& chain : model.chains) {
        if (chain_id != 0 && chain.name[0] != chain_id)
            continue;
        for (const auto& res : chain.residues) {
            for (const auto& atom : res.atoms) {
                float b = atom.b_iso;
                if (b >= 0.0f) {
                    // b_to_bin(b)
                    int bin_idx = static_cast<int>(
                        static_cast<double>(b) / static_cast<double>(hist.b_max) * hist.n_bins);
                    if (bin_idx >= hist.n_bins)
                        bin_idx = hist.n_bins - 1;
                    hist.b_vector[bin_idx].push_back(b);
                }
            }
        }
    }

    return hist;
}

} // namespace coot