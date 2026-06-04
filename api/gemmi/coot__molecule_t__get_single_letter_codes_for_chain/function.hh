#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

std::vector<std::pair<residue_spec_t, std::string>>
get_single_letter_codes_for_chain_gemmi(
  const gemmi::Model &model, const std::string &chain_id);

} // namespace coot