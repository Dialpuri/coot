#pragma once

#include <utility>
#include <gemmi/model.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__is_het_residue/gemmi/function.hh"

namespace coot {
namespace molecule_t {

std::pair<short int, int>
next_residue_number_in_chain_gemmi(const gemmi::Chain& chain,
                                   bool new_res_no_by_hundreds);

} // namespace molecule_t
} // namespace coot