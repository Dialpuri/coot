#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <utility>
#include <clipper/core/coords.h>

namespace coot {

struct graph_match_info_t {
  bool success = false;
  clipper::RTop_orth rtop;
  double dist_score = 1e20;
  std::vector<std::pair<std::pair<std::string, std::string>,
                        std::pair<std::string, std::string> > > matching_atom_names;
  int n_match = -99;
};

coot::graph_match_info_t graph_match_gemmi(
    const gemmi::Residue *res_moving,
    const gemmi::Residue *res_reference,
    bool apply_rtop_flag,
    bool match_hydrogens_also);

} // namespace coot