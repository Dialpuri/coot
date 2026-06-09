#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "lidia-core/lbg-molfile.hh"

namespace coot {

struct atom_selection_container_t_gemmi {
    gemmi::Structure* mol;
    int n_selected_atoms;
    int read_success;

    atom_selection_container_t_gemmi()
        : mol(nullptr), n_selected_atoms(0), read_success(0) {}
};

atom_selection_container_t_gemmi mdl_mol_to_asc_gemmi(
    const lig_build::molfile_molecule_t &m,
    float b_factor);

} // namespace coot