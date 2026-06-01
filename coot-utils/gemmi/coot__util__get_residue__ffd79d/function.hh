#pragma once

#include <string>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_residue_by_binary_search/gemmi/function.hh"

namespace coot {
namespace util {

inline gemmi::Residue *
get_residue_gemmi(const residue_spec_t &rs, gemmi::Structure &st) {
    return get_residue_by_binary_search_gemmi(rs.chain_id, rs.res_no, rs.ins_code, st);
}

} // namespace util
} // namespace coot