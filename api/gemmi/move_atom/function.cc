#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {

bool move_atom_gemmi(const std::string &atom_name_in, gemmi::Residue *res_p, const gemmi::Vec3 &new_pos) {
    // just change the position of the first atom that matches atom_name_in
    bool done = false;

    for (gemmi::Atom& atom : res_p->atoms) {
        if (atom.name == atom_name_in) {
            atom.pos.x = new_pos.x;
            atom.pos.y = new_pos.y;
            atom.pos.z = new_pos.z;
            done = true;
            break;
        }
    }
    return done;
}

} // namespace coot