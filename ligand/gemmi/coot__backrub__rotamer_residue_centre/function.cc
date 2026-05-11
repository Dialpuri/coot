#include "function.hh"
#include <gemmi/model.hpp>

namespace coot {
namespace backrub {

clipper::Coord_orth rotamer_residue_centre_gemmi(const gemmi::Residue& residue) {
    float sum_x = 0, sum_y = 0, sum_z = 0;
    int n_atoms = static_cast<int>(residue.atoms.size());
    
    for (const auto& atom : residue.atoms) {
        sum_x += atom.pos.x;
        sum_y += atom.pos.y;
        sum_z += atom.pos.z;
    }
    
    if (n_atoms > 0) {
        float inv = 1.0f / static_cast<float>(n_atoms);
        return clipper::Coord_orth(sum_x * inv, sum_y * inv, sum_z * inv);
    } else {
        return clipper::Coord_orth(0, 0, 0);
    }
}

} // namespace backrub
} // namespace coot