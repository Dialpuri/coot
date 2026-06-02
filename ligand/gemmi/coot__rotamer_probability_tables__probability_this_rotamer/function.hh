#pragma once

#include "ligand/rotamer.hh"
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>

namespace coot {

// Free-function port of the member function probability_this_rotamer.
// Takes the rotamer_probability_tables object as the first argument.
rotamer_probability_info_t probability_this_rotamer_gemmi(
  const rotamer_probability_tables &tables_obj,
  unsigned int i_table,
  const std::vector<std::pair<int,float>> &chi_angles);

} // namespace coot