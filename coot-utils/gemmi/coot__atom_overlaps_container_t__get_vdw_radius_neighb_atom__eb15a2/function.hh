#pragma once

#include <vector>

namespace coot {

class atom_overlaps_container_t {
public:
    std::vector<double> neighb_atom_radius;

    double get_vdw_radius_neighb_atom_gemmi(int idx_neigh_atom) const {
        // no index checking (a bit cowboy?)
        //
        double r = neighb_atom_radius[idx_neigh_atom];
        return r;
    }
};

} // namespace coot
