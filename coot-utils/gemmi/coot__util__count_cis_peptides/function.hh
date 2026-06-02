#pragma once

#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__cis_peptides_info_from_coords/gemmi/function.hh"

namespace coot {
namespace util {

inline int count_cis_peptides_gemmi(const gemmi::Structure &st) {
    return static_cast<int>(cis_peptides_info_from_coords_gemmi(st).size());
}

} // namespace util
} // namespace coot