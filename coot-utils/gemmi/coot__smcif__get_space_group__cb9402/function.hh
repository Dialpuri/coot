#pragma once
#include <string>
#include <vector>
#include <clipper/core/spacegroup.h>

namespace coot {
namespace smcif {

std::pair<bool, clipper::Spacegroup>
get_space_group_gemmi(const std::vector<std::string> &symm_strings);

} // namespace smcif
} // namespace coot