#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {
namespace util {

float average_temperature_factor_gemmi(
    std::vector<const gemmi::Atom*> atom_selection,
    int n_atoms,
    float low_cutoff,
    float high_cutoff,
    short int apply_low_cutoff,
    short int apply_high_cutoff) {

    float this_b = 0.0f;
    float b_sum = 0.0f;
    int n_sum = 0;

    for (int i = 0; i < n_atoms; i++) {
        this_b = atom_selection[i]->b_iso;
        if ((apply_low_cutoff && (this_b > low_cutoff)) || !apply_low_cutoff) {
            if ((apply_high_cutoff && (this_b > high_cutoff)) || !apply_high_cutoff) {
                b_sum += this_b;
                n_sum++;
            }
        }
    }

    float mean = 0.0f;
    if (n_atoms > 0)
        mean = b_sum / float(n_atoms);
    return mean;
}

} // namespace util
} // namespace coot
