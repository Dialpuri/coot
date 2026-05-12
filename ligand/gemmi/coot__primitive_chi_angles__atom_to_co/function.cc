#include "function.hh"

namespace coot {
namespace primitive_chi_angles {

clipper::Coord_orth atom_to_co_gemmi(const gemmi::Atom& atom) {
    return clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
}

} // namespace primitive_chi_angles
} // namespace coot