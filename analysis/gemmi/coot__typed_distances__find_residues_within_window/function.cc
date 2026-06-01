#include "function.hh"

namespace coot {
namespace typed_distances {

void find_residues_within_window_gemmi(
    const gemmi::Model& model, int half_wl,
    std::map<gemmi::CRA, std::vector<gemmi::CRA>, cra_less>& residues_within_window) {

  for (const gemmi::Chain& chain : model.chains) {
    int nres = static_cast<int>(chain.residues.size());
    for (int ires = 0; ires < nres; ires++) {
      const gemmi::Residue& residue = chain.residues[ires];
      gemmi::CRA res_cra{const_cast<gemmi::Chain*>(&chain), const_cast<gemmi::Residue*>(&residue), nullptr};
      for (int window_res = ires - half_wl; window_res <= ires + half_wl; window_res++) {
        if (window_res >= 0) {
          if (window_res < nres) {
            const gemmi::Residue& r = chain.residues[window_res];
            gemmi::CRA r_cra{const_cast<gemmi::Chain*>(&chain), const_cast<gemmi::Residue*>(&r), nullptr};
            residues_within_window[res_cra].push_back(r_cra);
          }
        }
      }
    }
  }
}

} // namespace typed_distances
} // namespace coot