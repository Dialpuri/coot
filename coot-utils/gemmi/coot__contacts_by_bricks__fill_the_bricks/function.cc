#include "function.hh"

namespace coot {
namespace contacts_by_bricks {

namespace {
inline unsigned int idx_3d_to_idx_1d(int idx_3d[3], int range[3]) {
    return (unsigned int)(
        idx_3d[0] +
        idx_3d[1] * range[0] +
        idx_3d[2] * range[0] * range[1]);
}
} // anonymous namespace

void fill_the_bricks_gemmi(
    gemmi::Atom const* atoms[],
    int n_atoms,
    float brick_size,
    const float lower_left[3],
    const int range[3],
    std::vector<std::set<unsigned int>>& atoms_in_bricks)
{
    float inv_brick_size = 1.0f / brick_size;
    for (int i = 0; i < n_atoms; i++) {
        gemmi::Atom const* at = atoms[i];
        int idx_3d[3];
        idx_3d[0] = static_cast<int>((at->pos.x - lower_left[0]) * inv_brick_size);
        idx_3d[1] = static_cast<int>((at->pos.y - lower_left[1]) * inv_brick_size);
        idx_3d[2] = static_cast<int>((at->pos.z - lower_left[2]) * inv_brick_size);
        unsigned int idx_1d = idx_3d_to_idx_1d(idx_3d, const_cast<int*>(range));
        // atoms that fly over the edge don't have NBCs :-)
        unsigned int n_bricks = atoms_in_bricks.size();
        if (idx_1d < n_bricks) {
            std::set<unsigned int>& ss = atoms_in_bricks.at(idx_1d);
            ss.insert(static_cast<unsigned int>(i));

            // now delete i from the neighbour brick sets
            // (disabled as in original)
            if (false) { // this block ~5ms for 400 residues
                for (int ix = idx_3d[0]-1; ix <= idx_3d[0]+1; ix++) {
                    if (ix >= 0) {
                        if (ix < range[0]) {
                            for (int iy = idx_3d[1]-1; iy <= idx_3d[1]+1; iy++) {
                                if (iy >= 0) {
                                    if (iy < range[1]) {
                                        for (int iz = idx_3d[2]-1; iz <= idx_3d[2]+1; iz++) {
                                            if (iz >= 0) {
                                                if (iz < range[2]) {
                                                    if (! ((ix==idx_3d[0]) && (iy==idx_3d[1]) && (iz==idx_3d[2]))) {
                                                        int idx_3d_neighb[3];
                                                        idx_3d_neighb[0] = ix;
                                                        idx_3d_neighb[1] = iy;
                                                        idx_3d_neighb[2] = iz;
                                                        unsigned int idx_neighb = idx_3d_to_idx_1d(idx_3d_neighb, const_cast<int*>(range));
                                                        if (idx_neighb < n_bricks)
                                                            atoms_in_bricks.at(idx_neighb).erase(static_cast<unsigned int>(i));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace contacts_by_bricks
} // namespace coot