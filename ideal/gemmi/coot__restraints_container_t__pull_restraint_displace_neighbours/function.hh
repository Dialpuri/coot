#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <set>
#include <cmath>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// gemmi port of coot::restraints_container_t::pull_restraint_displace_neighbours
// Displaces atoms within `radius` of `pull_atom` proportionally toward
// `new_pull_atom_target_position`, scaled by sqrt(1 - d/radius) and a factor of 0.5.
// Atoms at indices in `fixed_atom_indices` are skipped.
inline void pull_restraint_displace_neighbours_gemmi(
    gemmi::Atom* pull_atom,
    std::vector<gemmi::Atom*> &atom_ptr,
    const std::set<unsigned int> &fixed_atom_indices,
    const clipper::Coord_orth& new_pull_atom_target_position,
    float radius)
{
    bool use_top_hat_function = false;

    clipper::Coord_orth atom_current_position = co_gemmi(pull_atom);
    clipper::Coord_orth delta = new_pull_atom_target_position - atom_current_position;
    float f = 0.5f;
    float r_squared = radius * radius;

    for (unsigned int iat = 0; iat < atom_ptr.size(); iat++) {
        gemmi::Atom* at = atom_ptr[iat];
        if (fixed_atom_indices.find(iat) == fixed_atom_indices.end()) {
            float d_squared =
                (at->pos.x - pull_atom->pos.x) * (at->pos.x - pull_atom->pos.x) +
                (at->pos.y - pull_atom->pos.y) * (at->pos.y - pull_atom->pos.y) +
                (at->pos.z - pull_atom->pos.z) * (at->pos.z - pull_atom->pos.z);
            if (d_squared < r_squared) {
                float d = std::sqrt(d_squared);
                float ff = 1.0f - d / radius;
                if (ff < 0) ff = 0;
                float sf = std::sqrt(ff);
                if (use_top_hat_function)
                    sf = 1.0f;
                at->pos.x += sf * f * delta.x();
                at->pos.y += sf * f * delta.y();
                at->pos.z += sf * f * delta.z();
            }
        }
    }
}

} // namespace coot