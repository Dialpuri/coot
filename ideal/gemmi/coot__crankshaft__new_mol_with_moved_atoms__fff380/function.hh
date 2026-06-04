#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <vector>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__crankshaft_set__move_the_atoms/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__update_position/gemmi/function.hh"

namespace coot {

namespace crankshaft {

// Gemmi-port of scored_triple_angle_set_t (which is typedef'd to triple_crankshaft_set)
// Holds 3 crankshaft_sets and their corresponding rotation angles.
struct scored_triple_angle_set_t_gemmi {
    std::vector<coot::crankshaft_set> cs;
    std::vector<float> angles;

    scored_triple_angle_set_t_gemmi(
            gemmi::Residue& r0, gemmi::Residue& r1, gemmi::Residue& r2,
            gemmi::Residue& r3, gemmi::Residue& r4, gemmi::Residue& r5,
            const std::vector<float>& angles_in = {0.0f, 0.0f, 0.0f})
        : cs{coot::crankshaft_set(r0, r1, r2, r3),
             coot::crankshaft_set(r1, r2, r3, r4),
             coot::crankshaft_set(r2, r3, r4, r5)},
          angles(angles_in) {}

    coot::crankshaft_set& operator[](unsigned int i) { return cs[i]; }
    const coot::crankshaft_set& operator[](unsigned int i) const { return cs[i]; }
};

// Gemmi port of coot::crankshaft::new_mol_with_moved_atoms.
// Moves atoms in the crankshaft sets, copies the model (with moved atoms),
// then restores original positions in the source model. Returns the copy.
inline gemmi::Model new_mol_with_moved_atoms_gemmi(
    const gemmi::Model& mol,
    scored_triple_angle_set_t_gemmi& sas)
{
    std::map<gemmi::Atom*, clipper::Coord_orth> original_positions;

    int indices[] = { 2, 3, 4, 5 };
    for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t iat = 0; iat < 4; iat++) {
            gemmi::Atom* at = sas[i].v[indices[iat]];
            if (at) {
                clipper::Coord_orth pos = coot::co_gemmi(at);
                original_positions[at] = pos;
            }
        }
        sas[i].move_the_atoms_gemmi(sas.angles[i]);
    }

    // Copy the model (atoms now in moved positions)
    gemmi::Model mol_new = mol;

    // Restore original positions in the source model
    for (auto& kv : original_positions) {
        coot::update_position_gemmi(kv.first, kv.second);
    }

    return mol_new;
}

} // namespace crankshaft

} // namespace coot