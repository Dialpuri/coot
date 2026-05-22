#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <iostream>

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__get_first_residue/gemmi/function.hh"

namespace coot {
namespace ligand {

// Port of limit_solutions to gemmi
// This function limits solutions based on correlation threshold
void limit_solutions_gemmi(unsigned int iclust,
                          float frac_max_correl_lim,
                          int max_n_solutions,
                          float tolerance,
                          bool filter_by_torsion_match);

} // namespace ligand
} // namespace coot