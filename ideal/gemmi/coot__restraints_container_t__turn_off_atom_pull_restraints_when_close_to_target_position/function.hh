#pragma once

#include <vector>
#include <cmath>
#include <gemmi/model.hpp>
#include <coot/geometry/residue-and-atom-specs.hh>
#include <clipper/core/coords.h>

// Include for simple_restraint and restraint_type_t
// This is needed for the types used below
#include <coot/ideal/simple-restraint.hh>

namespace coot {

// gemmi port: free-function version of
// restraints_container_t::turn_off_atom_pull_restraints_when_close_to_target_position
//
// Scans restraints for TARGET_POS_RESTRAINT entries whose current atom position
// is within close_dist of the atom_pull_target_pos (excluding the dragged atom).
// Calls close() on matching restraints and returns their atom_specs.
inline std::vector<atom_spec_t> turn_off_atom_pull_restraints_when_close_to_target_position_gemmi(
    std::vector<simple_restraint>& restraints,
    const std::vector<const gemmi::Atom*>& atoms,
    const atom_spec_t& dragged_atom_spec)
{
    std::vector<atom_spec_t> v;
    double close_dist = 0.6;

    for (auto& r : restraints) {
        if (r.restraint_type == restraint_type_t(TARGET_POS_RESTRAINT)) {
            if (r.is_closed) {
                continue;
            } else {
                if (static_cast<size_t>(r.atom_index_1) >= atoms.size()) continue;
                const gemmi::Atom* at = atoms[r.atom_index_1];
                if (!at) continue;
                if (r.atom_spec != dragged_atom_spec) {
                    clipper::Coord_orth pos(at->pos.x, at->pos.y, at->pos.z);
                    double d = std::sqrt((pos - r.atom_pull_target_pos).lengthsq());
                    if (d < close_dist) {
                        r.close();
                        v.push_back(r.atom_spec);
                    }
                }
            }
        }
    }

    return v;
}

} // namespace coot