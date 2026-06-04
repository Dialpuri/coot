#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <string>

namespace coot {
namespace util {

// Returns {protein_count, nucleotide_count} for residues in the given chain.
// A null pointer is accepted and returns {0, 0}.
std::pair<unsigned int, unsigned int>
get_number_of_protein_or_nucleotides_gemmi(gemmi::Chain* chain_p);

} // namespace util
} // namespace coot