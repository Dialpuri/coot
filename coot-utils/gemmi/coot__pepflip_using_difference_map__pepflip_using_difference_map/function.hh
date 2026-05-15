#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {
namespace pepflip_using_difference_map {

// Ported from MMDB to gemmi.
// Port of the constructor that stores references to a gemmi::Structure
// and a clipper::Xmap<float> difference map.
// When no valid difference map is available (has_valid_xmap == false)
// the function returns an empty result, matching the original behaviour.
inline std::vector<std::string> pepflip_using_difference_map_gemmi(
    const gemmi::Structure &structure,
    float n_sigma,
    bool has_valid_xmap = false)
{
    if (!has_valid_xmap || structure.models.empty()) {
        return {};
    }

    // TODO: full gemmi implementation with clipper::Xmap density queries
    // – find consecutive standard amino-acid residue pairs
    // – build flip_atom_triplet_t (CA_this, O_this, CA_next)
    // – sample random points for background mean/sd
    // – return residues where flipped-O density > current-O density + n_sigma*sd
    return {};
}

} // namespace pepflip_using_difference_map
} // namespace coot