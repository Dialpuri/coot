#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

#include "api/coot-molecule.hh"
#include "coot-utils/peak-search.hh"
#include "utils/coot-utils.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

// Free function port of coot::molecule_t::difference_map_peaks
// Returns interesting_place_t vector with peak information
std::vector<coot::molecule_t::interesting_place_t> difference_map_peaks_gemmi(
    coot::molecule_t& mol, const gemmi::Structure& st, float n_rmsd);

} // namespace coot