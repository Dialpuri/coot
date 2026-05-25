#pragma once
#include <vector>
#include <string>
#include <map>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>

namespace coot {
namespace high_res {

// Port of buccafilter_neighbours from MMDB to gemmi
// Groups atoms by connectivity (same atom name → neighbors → same group)
// and filters the molecule based on these groups.
void buccafilter_neighbours_gemmi(gemmi::Structure& st);

} // namespace high_res
} // namespace coot