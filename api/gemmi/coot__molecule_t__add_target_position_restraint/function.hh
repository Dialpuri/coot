#pragma once

#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__molecule_t__cid_to_atom/gemmi/function.hh"

namespace coot {
namespace molecule_t {

// Add (or update) a target position restraint for an atom identified by its CID.
// If the atom already has a restraint, its target position is updated.
// Otherwise, a new entry is appended.
// Invalid CIDs are silently ignored.
void add_target_position_restraint_gemmi(
    const std::string &atom_cid,
    float pos_x, float pos_y, float pos_z,
    gemmi::Structure &st,
    std::vector<std::pair<gemmi::Atom *, clipper::Coord_orth>> &atoms_with_position_restraints);

} // namespace molecule_t
} // namespace coot