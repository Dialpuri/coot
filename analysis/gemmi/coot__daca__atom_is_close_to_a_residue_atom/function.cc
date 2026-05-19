#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {
namespace daca {

bool atom_is_close_to_a_residue_atom_gemmi(const gemmi::Atom& at, const gemmi::Residue& residue) {
    // Distance threshold: 1.7 + 1.7 + 1.5 = 4.9 Angstroms
    const float d_close = 1.7f + 1.7f + 1.5f;
    const float dd_close = d_close * d_close;
    
    for (const gemmi::Atom& ref_at : residue.atoms) {
        float dd =
            (at.pos.x - ref_at.pos.x) * (at.pos.x - ref_at.pos.x) +
            (at.pos.y - ref_at.pos.y) * (at.pos.y - ref_at.pos.y) +
            (at.pos.z - ref_at.pos.z) * (at.pos.z - ref_at.pos.z);
        if (dd < dd_close) {
            return true;
        }
    }
    return false;
}

}
}