#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/clipper_types.h"
#include "clipper/core/coords.h"

namespace coot {
namespace sequence_assignment {

class side_chain_score_t {
public:
  void move_std_res_to_this_res_pos_gemmi(const clipper::RTop_orth &rtop,
                                          gemmi::Residue *std_residue);
};

} // namespace sequence_assignment
} // namespace coot