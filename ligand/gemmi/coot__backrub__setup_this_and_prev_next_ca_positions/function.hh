#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include <string>
#include <stdexcept>

namespace coot {
namespace backrub {

void setup_this_and_prev_next_ca_positions_gemmi(
    const gemmi::Residue* this_res,
    const gemmi::Residue* prev_res,
    const gemmi::Residue* next_res,
    const std::string& alt_conf,
    clipper::Coord_orth& ca_this,
    clipper::Coord_orth& ca_prev,
    clipper::Coord_orth& ca_next
);

} // namespace backrub
} // namespace coot