#pragma once
#include <gemmi/model.hpp>

namespace coot {

// Deep copy a residue with its atoms.
// Takes CRA to allow reconstruction of chain context if needed.
// Returns a new Residue with copied atoms (chain is not set in returned residue).
gemmi::Residue deep_copy_residue_gemmi(gemmi::CRA cra);

}