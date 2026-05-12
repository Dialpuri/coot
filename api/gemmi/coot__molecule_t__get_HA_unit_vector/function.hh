#pragma once
#include <gemmi/model.hpp>
#include <utility>

namespace coot {

// Port of get_HA_unit_vector to gemmi
// Returns {true, unit_vector} if CA, C, N are found (CB optional)
// Returns {false, uninitialized Vec3} if required atoms are missing
std::pair<bool, gemmi::Vec3>
molecule_t_get_HA_unit_vector_gemmi(const gemmi::Residue& residue) {
    // Try to find atoms by name (gemmi uses clean names without padding)
    const gemmi::Atom* ca = nullptr;
    const gemmi::Atom* c = nullptr;
    const gemmi::Atom* n = nullptr;
    const gemmi::Atom* cb = nullptr;

    for (const gemmi::Atom& atom : residue.atoms) {
        if (atom.name == "CA") ca = &atom;
        else if (atom.name == "C") c = &atom;
        else if (atom.name == "N") n = &atom;
        else if (atom.name == "CB") cb = &atom;
    }

    if (ca && c && n) {
        gemmi::Vec3 ca_pos(ca->pos.x, ca->pos.y, ca->pos.z);
        gemmi::Vec3 c_pos(c->pos.x, c->pos.y, c->pos.z);
        gemmi::Vec3 n_pos(n->pos.x, n->pos.y, n->pos.z);
        
        gemmi::Vec3 dir_1 = ca_pos - c_pos;
        gemmi::Vec3 dir_2 = ca_pos - n_pos;
        
        gemmi::Vec3 r = dir_1 + dir_2;
        
        if (cb) {
            gemmi::Vec3 cb_pos(cb->pos.x, cb->pos.y, cb->pos.z);
            gemmi::Vec3 dir_3 = ca_pos - cb_pos;
            r = r + dir_3;
        }
        
        gemmi::Vec3 dir = r.normalized();
        return {true, dir};
    }
    
    return {false, gemmi::Vec3{0, 0, 0}};
}

} // namespace coot