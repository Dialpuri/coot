#pragma once

#include <gemmi/model.hpp>
#include <vector>

namespace coot {

std::vector<std::vector<const gemmi::Atom*>>
get_daca_fragments_gemmi(const gemmi::Residue* reference_residue_p);

} // namespace coot