#pragma once
#include <string>
#include <gemmi/model.hpp>

namespace coot {

// Ported from coot::get_term_type for gemmi API
// Takes CRA (Chain-Residue-Atom) to access parent context since gemmi Residue/Atom have no parent pointer
std::string get_term_type_gemmi(const gemmi::CRA& cra);

} // namespace coot