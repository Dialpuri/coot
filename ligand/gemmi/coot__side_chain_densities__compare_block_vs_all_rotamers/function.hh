#pragma once

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/core/xmap.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__side_chain_densities__compare_block_vs_rotamer/gemmi/function.hh"
#include "utils/coot-utils.hh"

namespace coot {

namespace side_chain_densities {

std::map<std::string, std::pair<std::string, double>>
compare_block_vs_all_rotamers_gemmi(
    density_box_t block,
    gemmi::Residue *residue_p,
    const std::string &data_dir,
    const std::pair<bool, std::vector<std::pair<std::string, std::string>>> &rotamer_limits,
    const clipper::Xmap<float> &xmap);

} // namespace side_chain_densities

} // namespace coot
