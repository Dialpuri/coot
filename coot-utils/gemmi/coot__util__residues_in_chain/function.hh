#pragma once

#include <vector>
#include <string>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Return all residues in a chain with the given chain ID
std::vector<gemmi::Residue*> residues_in_chain_gemmi(gemmi::Model& model, const std::string& chain_id);

}} // namespace coot::util