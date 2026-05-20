#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

bool residue_sorter_gemmi(const std::pair<bool, gemmi::CRA> &r1,
                          const std::pair<bool, gemmi::CRA> &r2);

} // namespace coot