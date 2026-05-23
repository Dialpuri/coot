#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <coot/ligand/rotamer.hh>

namespace coot { namespace rotamer_helper {

// Ported from MMDB version that took mmdb::Residue*
// Uses gemmi::CRA to access parent context (chain, residue, atom)
rotamer_probability_info_t probability_of_this_rotamer_gemmi(const gemmi::CRA& cra);

} // namespace rotamer_helper
} // namespace coot