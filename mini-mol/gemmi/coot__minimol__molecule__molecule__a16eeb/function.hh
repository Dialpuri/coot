#pragma once

#include <gemmi/model.hpp>
#include "mini-mol/mini-mol.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__molecule__setup/gemmi/function.hh"

namespace coot {
namespace minimol {

// gemmi port of the molecule(mmdb::Manager*, bool) constructor.
// Creates a molecule from a gemmi::Structure (factory function since
// we cannot add constructors to the existing molecule class).
inline molecule molecule_gemmi(const gemmi::Structure* st, bool udd_atom_index_to_user_data_flag = false) {
    molecule mol;
    setup_gemmi(mol, st);
    // udd_atom_index_to_user_data_flag is MMDB-only; ignored in gemmi
    (void)udd_atom_index_to_user_data_flag;
    return mol;
}

}  // namespace minimol
}  // namespace coot