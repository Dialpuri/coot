#pragma once
#include <utility>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__cho__is_het_residue/gemmi/function.hh"

namespace coot {
namespace cho {

// Port of coot::cho::next_residue_number_in_chain from MMDB to gemmi
// Returns {state, next_residue_number} where state=1 means success, state=0 means no valid chain
std::pair<short int, int> next_residue_number_in_chain_gemmi(const gemmi::Chain* chain, bool new_res_no_by_hundreds);

} // namespace cho
} // namespace coot