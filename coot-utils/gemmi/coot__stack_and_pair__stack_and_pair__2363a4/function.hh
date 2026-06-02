#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <map>
#include <vector>
#include <utility>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__stack_and_pair__calculate_residue_normals__bc3097/gemmi/function.hh"

namespace coot {
namespace stack_and_pair {

/// Gemmi-port of the stack_and_pair constructor:
///  init() + calculate_residue_normals(residues_vec)
///
/// Returns a struct holding the populated normal_map.
struct stack_and_pair_result {
    std::map<const gemmi::Residue*, clipper::Coord_orth> normal_map;
};

inline stack_and_pair_result stack_and_pair_gemmi(const std::vector<std::pair<bool, const gemmi::Residue*>>& residues_vec) {
    stack_and_pair_result result;
    result.normal_map = calculate_residue_normals_gemmi(residues_vec);
    return result;
}

} // namespace stack_and_pair
} // namespace coot