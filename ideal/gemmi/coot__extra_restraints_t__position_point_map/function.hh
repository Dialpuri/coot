#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <map>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

class extra_restraints_t {
public:
    std::map<gemmi::Atom*, clipper::Coord_orth> position_point_map_gemmi(
        const gemmi::Structure* st_running,
        const gemmi::Structure* st_ref) const;
};

} // namespace coot
