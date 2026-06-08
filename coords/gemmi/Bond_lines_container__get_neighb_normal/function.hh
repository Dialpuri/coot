#pragma once

#include "clipper/core/coords.h"
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "coot/coot-utils/coot-coord-utils.hh"

#include <vector>
#include <string>
#include <algorithm>

/// gemmi port of Bond_lines_container::get_neighb_normal
/// Computes the normal to the plane defined by neighbours of two bonded atoms.
clipper::Coord_orth get_neighb_normal_gemmi(
    const coot::protein_geometry& geom,
    int imol,
    int iat_1,
    int iat_2,
    const std::vector<gemmi::Atom>& atoms,
    const std::string& res_name,
    bool also_2nd_order_neighbs_flag)
{
    clipper::Coord_orth pt(0, 0, 0);

    std::string at_n_1 = coot::atom_id_mmdb_expand(atoms[iat_1].name);
    std::string at_n_2 = coot::atom_id_mmdb_expand(atoms[iat_2].name);

    std::vector<std::string> neighbours = geom.get_bonded_neighbours(
        res_name, imol, at_n_1, at_n_2, also_2nd_order_neighbs_flag);

    char alt_conf_bond = atoms[iat_1].altloc;

    if (neighbours.size() > 2) {
        // Collect neighbour atom indices (avoid duplicates by index)
        std::vector<int> neighb_atom_indices;
        for (unsigned int i = 0; i < neighbours.size(); i++) {
            for (unsigned int j = 0; j < atoms.size(); j++) {
                std::string atom_name = coot::atom_id_mmdb_expand(atoms[j].name);
                if (neighbours[i] == atom_name) {
                    if (atoms[j].altloc == alt_conf_bond) {
                        // only add them if they are not there already (belt and braces test)
                        if (std::find(neighb_atom_indices.begin(),
                                      neighb_atom_indices.end(), j) ==
                            neighb_atom_indices.end()) {
                            neighb_atom_indices.push_back(static_cast<int>(j));
                        }
                    }
                }
            }
        }

        if (neighb_atom_indices.size() > 2) {
            std::vector<clipper::Coord_orth> neighb_atoms_pos(neighb_atom_indices.size());
            for (unsigned int i = 0; i < neighb_atom_indices.size(); i++) {
                int idx = neighb_atom_indices[i];
                neighb_atoms_pos[i] = clipper::Coord_orth(
                    atoms[idx].pos.x,
                    atoms[idx].pos.y,
                    atoms[idx].pos.z);
            }
            coot::lsq_plane_info_t lp(neighb_atoms_pos);
            pt = lp.normal();
        }
    }

    if (pt.x() == 0 && pt.y() == 0 && pt.z() == 0) {
        // make something up...
        pt = clipper::Coord_orth(0, 0, 1);
    }

    return pt;
}
