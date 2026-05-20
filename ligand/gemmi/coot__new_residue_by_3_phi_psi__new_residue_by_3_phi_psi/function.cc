#include "function.hh"
#include <gemmi/pdb.hpp>

namespace coot {

new_residue_by_3_phi_psi::new_residue_by_3_phi_psi(
    const std::string &terminus_type_in,
    const gemmi::CRA &cra,
    const std::vector<gemmi::Residue> &all_residues)
    : cra(cra),
      terminus_type(terminus_type_in),
      n_threads(0),
      thread_pool_p(nullptr),
      upstream_neighbour_residue_p(nullptr),
      downstream_neighbour_residue_p(nullptr) {

    chain_id = cra.chain->name;
    init_phi_psi_plot();
}

void new_residue_by_3_phi_psi::init_phi_psi_plot() {
    upstream_neighbour_residue_p = nullptr;
    downstream_neighbour_residue_p = nullptr;
}

} // namespace coot