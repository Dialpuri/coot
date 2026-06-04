#pragma once

#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

class pepflip_using_difference_map {
public:
    pepflip_using_difference_map(gemmi::Model* model);
    std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth>>
    make_random_other_pairs_gemmi(int n_others) const;

private:
    gemmi::Model* model_;
};

} // namespace coot