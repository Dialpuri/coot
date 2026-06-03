#pragma once
#include "coot-utils/tree.hh"
#include <gemmi/model.hpp>
#include <string>

namespace coot {

struct linked_residue_t_gemmi {
  gemmi::Residue* residue;
  std::string link_type;
};

// gemmi port of glyco_tree_t::get_level — free function because C++ does not
// allow re-opening coot::glyco_tree_t (which is mmdb-based).
int get_level_gemmi(const tree<linked_residue_t_gemmi>& glyco_tree,
                    gemmi::Residue* residue_p);

} // namespace coot