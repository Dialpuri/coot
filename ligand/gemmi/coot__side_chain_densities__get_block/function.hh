#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <iostream>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__density_box_t__density_box_t/gemmi/function.hh"

namespace coot {

namespace side_chain_densities {

inline density_box_t get_block_gemmi(
    const std::map<gemmi::Residue*, density_box_t>& density_block_map_cache,
    gemmi::Residue* residue_p) {

    auto it = density_block_map_cache.find(residue_p);
    if (it == density_block_map_cache.end()) {
        std::cout << "ERROR:: in get_block(): Hideous failure!" << std::endl;
    }
    return it->second;
}

} // namespace side_chain_densities
} // namespace coot