#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "lidia-core/lbg-shared.hh"

namespace pli {

static inline double get_radius_for_gemmi(const std::string &ele) {
    double radius = 1.70;
    if (ele == "H")  radius = 1.20;
    if (ele == "N")  radius = 1.55;
    if (ele == "O")  radius = 1.52;
    if (ele == "S")  radius = 1.80;
    return radius;
}

inline coot::bash_distance_t find_bash_distance_gemmi(
    const clipper::Coord_orth &ligand_atom_pos,
    const clipper::Coord_orth &hydrogen_pos,
    const std::vector<gemmi::CRA> &close_residue_atoms)
{
    double cannonball_radius = 0.8;
    double max_dist = 4.05;

    clipper::Coord_orth h_vector((hydrogen_pos - ligand_atom_pos).unit());

    // set the atomic radii
    std::vector<double> radius(close_residue_atoms.size());
    for (size_t iat = 0; iat < close_residue_atoms.size(); iat++) {
        const gemmi::Atom *atom = close_residue_atoms[iat].atom;
        std::string ele = atom ? atom->element.name() : "";
        radius[iat] = get_radius_for_gemmi(ele);
    }

    coot::bash_distance_t dd; // default: dist=-1, limited=false

    // cache atom positions
    std::vector<clipper::Coord_orth> atom_positions(close_residue_atoms.size());
    for (size_t i = 0; i < close_residue_atoms.size(); i++) {
        const gemmi::Atom *atom = close_residue_atoms[i].atom;
        if (atom)
            atom_positions[i] = clipper::Coord_orth(atom->pos.x, atom->pos.y, atom->pos.z);
    }

    for (double slide = 0; slide <= max_dist; slide += 0.04) {
        clipper::Coord_orth test_pt = ligand_atom_pos + slide * h_vector;
        for (size_t iat = 0; iat < atom_positions.size(); iat++) {
            double atom_radius_plus_cbr = radius[iat] + cannonball_radius;
            double d_squared = (test_pt - atom_positions[iat]).lengthsq();
            if (d_squared < atom_radius_plus_cbr * atom_radius_plus_cbr) {
                dd = coot::bash_distance_t(slide);
                break;
            }
        }
        if (dd.limited) {
            break;
        }
    }
    return dd;
}

} // namespace pli