#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class new_residue_by_3_phi_psi {
public:
    new_residue_by_3_phi_psi(const std::string &terminus_type_in,
                                   const gemmi::CRA &cra,
                                   const std::vector<gemmi::Residue> &all_residues);

    std::string chain_id;
    std::string terminus_type;
    gemmi::CRA cra;
    int n_threads;
    void *thread_pool_p;
    const gemmi::Residue *upstream_neighbour_residue_p;
    const gemmi::Residue *downstream_neighbour_residue_p;

private:
    void init_phi_psi_plot();
};

} // namespace coot