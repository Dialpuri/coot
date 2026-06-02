#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include <utility>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>

namespace coot {
namespace hole {

// make_atom_selection_gemmi: select atoms within 'radius' of pt,
// excluding HOH and keeping only H,C,O,N,S elements.
inline std::vector<gemmi::CRA>
make_atom_selection_gemmi(const gemmi::Model& model,
                          const clipper::Coord_orth& pt,
                          double radius) {
    std::vector<gemmi::CRA> result;
    double r_sq = radius * radius;
    for (const auto& chain : model.chains) {
        for (const auto& res : chain.residues) {
            if (res.name == "HOH") continue;
            for (const auto& atom : res.atoms) {
                const auto& elem = atom.element;
                if (!(elem.name() == "H" || elem.name() == "C" ||
                      elem.name() == "O" || elem.name() == "N" ||
                      elem.name() == "S")) continue;
                double dx = atom.pos.x - pt.x();
                double dy = atom.pos.y - pt.y();
                double dz = atom.pos.z - pt.z();
                if (dx*dx + dy*dy + dz*dz <= r_sq) {
                    result.push_back({const_cast<gemmi::Chain*>(&chain),
                                      const_cast<gemmi::Residue*>(&res),
                                      const_cast<gemmi::Atom*>(&atom)});
                }
            }
        }
    }
    return result;
}

// sphere_size_gemmi: minimum distance from pt to atoms in the selection.
// Returns -1.0 if no atoms are selected.
inline double sphere_size_gemmi(const clipper::Coord_orth& pt,
                                const std::vector<gemmi::CRA>& atom_list) {
    if (atom_list.empty()) return -1.0;
    double min_dist = std::numeric_limits<double>::max();
    for (const auto& cra : atom_list) {
        if (!cra.atom) continue;
        double dx = cra.atom->pos.x - pt.x();
        double dy = cra.atom->pos.y - pt.y();
        double dz = cra.atom->pos.z - pt.z();
        double d = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (d < min_dist) min_dist = d;
    }
    return min_dist;
}

// optimize_point_gemmi: gemmi port of coot::hole::optimize_point
//
// By moving about in the plane perpendicular to v_hat we try to maximise
// the distance between the point and the atoms around it.
// Returns the optimised position and the max distance (sphere_size).
inline std::pair<clipper::Coord_orth, double>
optimize_point_gemmi(const clipper::Coord_orth& pt,
                     const clipper::Coord_orth& v_hat,
                     const gemmi::Model& model,
                     double radius) {

    int max_moves_without_improvement = 80;
    int move_count = 0;
    double biggest_allowed_distance_to_protein = 5.0;

    // Build atom selection using gemmi equivalent
    std::vector<gemmi::CRA> atom_list =
        make_atom_selection_gemmi(model, pt, radius);

    clipper::Coord_orth current_pt = pt;
    double D_max = 0.1;
    float rmi = 1.0f / static_cast<float>(RAND_MAX);

    double current_ss = sphere_size_gemmi(pt, atom_list);

    while (move_count < max_moves_without_improvement) {
        double d1 = 2.0 * std::rand() * rmi - 1.0;
        double d2 = 2.0 * std::rand() * rmi - 1.0;
        double d3 = 2.0 * std::rand() * rmi - 1.0;
        clipper::Coord_orth y_rand(d1, d2, d3);
        clipper::Coord_orth y_hat_rand(y_rand.unit());
        clipper::Coord_orth y_prime_rand =
            y_hat_rand - (clipper::Coord_orth::dot(v_hat, y_hat_rand)) * v_hat;

        clipper::Coord_orth trial_pt = current_pt + D_max * y_prime_rand;

        double ss = sphere_size_gemmi(trial_pt, atom_list);
        if (ss > current_ss) {
            current_pt = trial_pt;
            current_ss = ss;
            move_count = 0; // reset
            if (ss > biggest_allowed_distance_to_protein) {
                break; // the point has escaped
            }
        } else {
            move_count++;
        }
    }
    return std::make_pair(current_pt, current_ss);
}

} // namespace hole
} // namespace coot