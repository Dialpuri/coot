#pragma once
#include "coot-utils/glyco-tree.hh"
#include "coot/geometry/residue-and-atom-specs.hh"
#include <gemmi/model.hpp>
#include <fstream>
#include <string>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// gemmi port of glyco_tree_t::output_internal_distances
void output_internal_distances_gemmi(
    const gemmi::Residue &residue,
    const gemmi::Residue *parent_p,
    double dist_crit,
    std::ofstream &f,
    int model_number,
    const std::string &chain_id);

} // namespace coot