#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

#include "clipper/core/coords.h"

namespace coot {

std::string convert_base_name(const std::string &std_base_name, bool use_old_style_naming);

}

namespace coot {
namespace util {

void mutate_base_gemmi(
  gemmi::Residue &residue,
  const gemmi::Residue &std_base,
  bool use_old_style_naming,
  bool print_match_stats_flag,
  float b_factor);

}
}