#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot::util::map_to_model_correlation_stats_per_residue_run {

struct residue_run_t {
  std::vector<gemmi::Residue*> residues;
  void add_residue_gemmi(gemmi::Residue* r) { residues.push_back(r); }
};

} // namespace coot::util::map_to_model_correlation_stats_per_residue_run