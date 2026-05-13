#pragma once
#include <vector>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {

// Port of coot::co from MMDB to gemmi
// Takes a gemmi::Atom pointer and returns clipper::Coord_orth with its coordinates
inline clipper::Coord_orth co_gemmi(const gemmi::Atom* at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

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

    bool is_inside_an_env_atom_to_which_its_bonded_gemmi(
        int idx,
        const std::vector<int> &bonded_neighb_indices,
        const std::vector<gemmi::Atom*> &env_residue_atoms,
        const clipper::Coord_orth &pt_at_surface) const {
        bool r = false;
        double r_1 = get_vdw_radius_neighb_atom_gemmi(idx);
        for (size_t i = 0; i < bonded_neighb_indices.size(); i++) {
            gemmi::Atom *env_atom = env_residue_atoms[bonded_neighb_indices[i]];
            clipper::Coord_orth pt_env_atom = co_gemmi(env_atom);
            double r_2 = 1.6;
            std::string ele = env_atom->element.name();
            if (ele == "H") {
                r_2 = 0.97;
            }
            double r_2_sqrd = r_2 * r_2;
            clipper::Coord_orth diff = pt_at_surface - pt_env_atom;
            double d_sqrd = diff.lengthsq();
            if (d_sqrd < r_2_sqrd) {
                r = true;
                break;
            }
        }
        return r;
    }
};

} // namespace coot