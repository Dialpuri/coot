#pragma once

#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

/// Port of coot::get_selection_handle to gemmi.
/// Simulates MMDB selection handle allocation: returns -1 for nullptr,
/// otherwise allocates a handle starting from 2.
int get_selection_handle_gemmi(gemmi::Structure* mol, const coot::atom_spec_t& at);

/// Simulate mol->DeleteSelection() — frees the handle so it can be reused.
void delete_selection_gemmi(int sel_hnd);

} // namespace coot