#pragma once

#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <string>
#include <utility>
#include "coot/geometry/residue-and-atom-specs.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__link_atoms__96989c/gemmi/function.hh"

namespace coot {

class atom_overlaps_container_t {
public:
  bool is_linked_gemmi(gemmi::CRA cra1, gemmi::CRA cra2,
                       gemmi::Model* model,
                       gemmi::Structure* structure) const;
};

} // namespace coot