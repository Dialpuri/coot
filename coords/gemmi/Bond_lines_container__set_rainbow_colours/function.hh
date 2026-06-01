#pragma once

#include <gemmi/model.hpp>
#include "utils/coot-utils.hh"
#include <vector>

// Gemmi port of Bond_lines_container::set_rainbow_colours.
// Returns the UDD handle the MMDB version would register (16777217 on first call).
// Also collects rainbow values for each qualifying atom.
std::vector<float> set_rainbow_colours_gemmi(const gemmi::Structure& st);