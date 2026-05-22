#pragma once

#include <set>
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

// Port of coot::interface_residues to gemmi.
// Takes a gemmi Structure, two chain IDs, and a minimum distance threshold.
// Returns a pair of sets: residues in chain_A that interface with chain_B,
// and residues in chain_B that interface with chain_A.
std::pair<std::set<gemmi::Residue*>, std::set<gemmi::Residue*>>
interface_residues_gemmi(gemmi::Structure& st,
                          const std::string& chain_A,
                          const std::string& chain_B,
                          float min_dist);

} // namespace coot