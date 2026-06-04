#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <map>
#include <vector>
#include <utility>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

/// gemmi-port of atom_overlaps_container_t — only the members needed by
/// is_inside_another_ligand_atom_gemmi are declared here.
struct atom_overlaps_container_t {
    std::map<int, std::vector<std::pair<const gemmi::Atom*, double>>> ligand_atom_neighbour_map;
    double probe_radius = 1.4;

    /// gemmi port of is_inside_another_ligand_atom.
    /// Checks whether dot_pt lies inside the expanded vdw+probe sphere of
    /// any neighbour atom stored under the given idx.
    bool is_inside_another_ligand_atom_gemmi(int idx,
                                             const clipper::Coord_orth &probe_pos,
                                             const clipper::Coord_orth &dot_pt) const {
        bool r = false;

        if (idx >= 0) {
            auto it = ligand_atom_neighbour_map.find(idx);
            if (it != ligand_atom_neighbour_map.end()) {
                const auto &v = it->second;
                for (unsigned int i = 0; i < v.size(); i++) {
                    clipper::Coord_orth pt = coot::co_gemmi(v[i].first);
                    double dist_sqrd = (dot_pt - pt).lengthsq();
                    double radius_other = v[i].second;
                    radius_other += probe_radius;
                    if (dist_sqrd < radius_other * radius_other) {
                        r = true;
                        break;
                    }
                }
            }
        }
        return r;
    }
};

} // namespace coot