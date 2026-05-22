#pragma once
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Returns: -1 if residue is empty, 0 if no hetatms, 1 if hetatms present
int residue_has_hetatms_gemmi(const gemmi::Residue &res);

}} // namespace coot::util