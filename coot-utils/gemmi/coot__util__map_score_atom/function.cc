#include "function.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {
namespace util {

float map_score_atom_gemmi(const gemmi::Atom *atom,
                           const clipper::Xmap<float> &xmap) {
    float f = 0;
    if (atom) {
        f = density_at_point(xmap,
                             clipper::Coord_orth(atom->pos.x, atom->pos.y, atom->pos.z));
    }
    return f;
}

} // namespace util
} // namespace coot