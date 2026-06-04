#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <algorithm>

namespace coot {
namespace util {

inline float median_temperature_factor_gemmi(
    std::vector<gemmi::Atom const*> atom_selection,
    float low_cutoff,
    float high_cutoff,
    bool apply_low_cutoff,
    bool apply_high_cutoff) {

  float median = 0;
  std::vector<float> b;

  for (auto const* atom : atom_selection) {
    float this_b = atom->b_iso;
    if ((apply_low_cutoff && (this_b > low_cutoff)) || !apply_low_cutoff) {
      if ((apply_high_cutoff && (this_b > high_cutoff)) || !apply_high_cutoff) {
        b.push_back(this_b);
      }
    }
  }

  if (!b.empty()) {
    std::sort(b.begin(), b.end());
    int mid_pos = static_cast<int>(b.size()) / 2;
    median = b[mid_pos];
  }
  return median;
}

} // namespace util
} // namespace coot