#pragma once

#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include <utility>
#include "api/rama-plot-phi-psi.hh"

namespace rama_plot {
namespace util {

std::pair<bool, rama_plot::phi_psi_t>
get_phi_psi_gemmi(gemmi::CRA residue_0, gemmi::CRA residue_1, gemmi::CRA residue_2);

} // namespace util
} // namespace rama_plot