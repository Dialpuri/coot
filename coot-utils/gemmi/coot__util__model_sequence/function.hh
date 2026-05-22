#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>

namespace coot { namespace util {

// Convert a vector of gemmi residues to a sequence string.
// Mimics the original MMDB version: returns 1-letter codes for standard
// amino acids (or all residues if allow_ligands is true).
std::string model_sequence_gemmi(const std::vector<std::pair<gemmi::Residue*, int>>& sa,
                                  bool allow_ligands);

// Helper: check if residue name is a standard amino acid (including MSE)
bool is_standard_amino_acid_name(const std::string& residue_name);

}} // namespace coot::util