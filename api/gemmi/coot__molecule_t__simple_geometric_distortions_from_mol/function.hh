#pragma once

#include <string>
#include <utility>
#include <gemmi/model.hpp>

namespace coot {

// Port of molecule_t::simple_geometric_distortions_from_mol to gemmi
// Takes a gemmi::CRA instead of mmdb::Residue*
// Returns {status, distortion} pair
std::pair<int, double>
simple_geometric_distortions_from_mol_gemmi(const gemmi::CRA &cra,
                                           bool with_nbcs);

}