#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <cmath>
#include <string>
#include <vector>

#include "geometry/residue-and-atom-specs.hh"
#include "api/coot-molecule.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_residue/gemmi/function.hh"

namespace coot {

std::vector<atom_distance_t> get_distances_between_atoms_of_residues_gemmi(
    const std::string &cid_res_1,
    const std::string &cid_res_2,
    float dist_max,
    gemmi::Model &model);

} // namespace coot