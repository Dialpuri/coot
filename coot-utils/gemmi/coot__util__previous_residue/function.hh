#pragma once

#include <gemmi/model.hpp>
#include "coot-utils/coot-coord-utils.hh"

namespace coot {
namespace util {

gemmi::Residue *
previous_residue_gemmi(const gemmi::CRA& cra);

} // namespace util
} // namespace coot