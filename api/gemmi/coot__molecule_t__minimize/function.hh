#pragma once
#include <string>

namespace coot {
  class molecule_t {
  public:
    int minimize_gemmi(const std::string &atom_selection_cid,
                       int n_cycles,
                       bool do_rama_plot_restraints, float rama_plot_weight,
                       bool do_torsion_restraints, float torsion_weight, bool refinement_is_quiet,
                       void *geom_p);
  };
}