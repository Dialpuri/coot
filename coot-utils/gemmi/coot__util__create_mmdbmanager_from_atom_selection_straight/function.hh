#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>
#include <functional>

namespace coot {
namespace util {

// Ported from MMDB version: creates a new gemmi::Structure from atom selection
// The selection_callback is called for each atom in the input structure
// and should return true to include the atom in the output
gemmi::Structure create_mmdbmanager_from_atom_selection_straight_gemmi(
    const gemmi::Structure& orig_st,
    const std::function<bool(const gemmi::CRA&)>& selection_callback);

} // namespace util
} // namespace coot