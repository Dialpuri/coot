#pragma once
#include <gemmi/model.hpp>
#include <cmath>

namespace coot {
namespace util {

// Port of MMDB's standard_deviation_temperature_factor to gemmi
float standard_deviation_temperature_factor_gemmi(
    const gemmi::Model& model,
    float low_cutoff = 0.0f,
    float high_cutoff = 999.0f,
    bool apply_low_cutoff = false,
    bool apply_high_cutoff = false) {

    double b_sum = 0.0;
    double b_sum_sqs = 0.0;
    int n_sum = 0;
    int total_atoms = 0;

    for (auto& chain : model.chains) {
        for (auto& residue : chain.residues) {
            for (const auto& atom : residue.atoms) {
                total_atoms++;
                float this_b = atom.b_iso;
                if ((apply_low_cutoff && (this_b > low_cutoff)) ||
                    !apply_low_cutoff) {
                    if ((apply_high_cutoff && (this_b > high_cutoff)) ||
                        !apply_high_cutoff) {
                        b_sum += this_b;
                        b_sum_sqs += this_b * this_b;
                        n_sum++;
                    }
                }
            }
        }
    }

    double mean = 0.0;
    double var = 0.0;
    float sd = 0.0;
    if (total_atoms > 0) {
        mean = b_sum / double(total_atoms);
        var = b_sum_sqs / double(total_atoms) - mean * mean;
        if (var < 0.0)
            var = 0.0;
        sd = sqrt(var);
    }
    return sd;
}

} // namespace util
} // namespace coot