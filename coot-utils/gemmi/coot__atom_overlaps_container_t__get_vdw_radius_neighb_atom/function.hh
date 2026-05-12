#pragma once
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class atom_overlaps_container_t {
public:
    std::vector<double> neighb_atom_radius;
    
    double get_vdw_radius_neighb_atom_gemmi(int idx_neigh_atom) const {
        if (idx_neigh_atom < 0) {
            return 0.0;
        }
        if (static_cast<size_t>(idx_neigh_atom) >= neighb_atom_radius.size()) {
            if (neighb_atom_radius.empty()) {
                return 0.0;
            }
            return neighb_atom_radius.back();
        }
        return neighb_atom_radius[idx_neigh_atom];
    }
};

} // namespace coot