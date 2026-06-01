#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>

namespace coot {
namespace typed_distances {

// Comparator for gemmi::CRA to use as std::map key.
// Compares by (chain, residue, atom) pointer ordering.
struct cra_less {
  bool operator()(const gemmi::CRA& a, const gemmi::CRA& b) const {
    if (a.chain != b.chain) return a.chain < b.chain;
    if (a.residue != b.residue) return a.residue < b.residue;
    return a.atom < b.atom;
  }
};

// For each residue in the model, find all residues within a window of
// half_wl positions (by index) in the same chain.
// Appends to the given map (does NOT clear it first — matches original MMDB behavior).
void find_residues_within_window_gemmi(
    const gemmi::Model& model, int half_wl,
    std::map<gemmi::CRA, std::vector<gemmi::CRA>, cra_less>& residues_within_window);

} // namespace typed_distances
} // namespace coot