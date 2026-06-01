#pragma once

#include <gemmi/model.hpp>
#include "coot-utils/fragment-container.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

fragment_container_t make_overlapping_fragments_gemmi(
    const gemmi::Structure &st,
    const std::string &chain_id,
    unsigned int fragment_length);

} // namespace coot