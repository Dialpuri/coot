#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <utility>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__distance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__angle/gemmi/function.hh"

namespace coot {

// Gemmi-compatible h_bond struct
struct h_bond_gemmi {
    double dist;
    double angle_1;
    double angle_2;
    double angle_3;

    const gemmi::Atom *hb_hydrogen;
    const gemmi::Atom *acceptor;
    const gemmi::Atom *donor;
    const gemmi::Atom *donor_neigh;
    const gemmi::Atom *acceptor_neigh;

    int ligand_atom_is_donor;
    int hb_type;
    int hb_type_udd;

    h_bond_gemmi(const gemmi::Atom *at_h, const gemmi::Atom *at_a, int is_donor = 0)
        : dist(0.), angle_1(0.), angle_2(-1.), angle_3(-1.),
          hb_hydrogen(at_h), acceptor(at_a), donor(nullptr),
          donor_neigh(nullptr), acceptor_neigh(nullptr),
          ligand_atom_is_donor(is_donor), hb_type(0), hb_type_udd(-1) {}
};

namespace h_bonds {

std::pair<bool, h_bond_gemmi>
make_h_bond_from_ligand_hydrogen_gemmi(
    const gemmi::Atom *at_1,  // H on ligand
    const gemmi::Atom *at_2,  // acceptor on residue
    const std::vector<std::pair<const gemmi::Atom *, float>> &nb_1,
    const std::vector<std::pair<const gemmi::Atom *, float>> &nb_2) {

    h_bond_gemmi bond(at_1, at_2, 1); // ligand atom is Hydrogen
    bond.dist = coot::distance_gemmi(at_1, at_2);
    bool neighbour_distances_and_angles_are_good = true;
    bool good_donor_acceptor_dist = false;

    // Angle D-H-A
    for (unsigned int iD = 0; iD < nb_1.size(); iD++) {
        double angle = coot::angle_gemmi(nb_1[iD].first, at_1, at_2);
        double dist  = coot::distance_gemmi(nb_1[iD].first, at_2);
        if (dist < 3.9)
            good_donor_acceptor_dist = true;
        if (!bond.donor) {
            bond.donor = nb_1[iD].first;
            bond.angle_1 = angle;
        }
        if (angle < 90) {
            neighbour_distances_and_angles_are_good = false;
            break;
        }
    }

    // Angle H-A-AA
    for (unsigned int iA = 0; iA < nb_2.size(); iA++) {
        double angle = coot::angle_gemmi(at_1, at_2, nb_2[iA].first);
        if (!bond.acceptor) {
            bond.angle_2 = angle;
        }
        if (angle < 90) {
            neighbour_distances_and_angles_are_good = false;
            break;
        }
    }

    // Angle D-A-AA
    for (unsigned int iD = 0; iD < nb_1.size(); iD++) {
        for (unsigned int iA = 0; iA < nb_2.size(); iA++) {
            double angle = coot::angle_gemmi(nb_1[iD].first, at_2, nb_2[iA].first);
            if (!bond.acceptor_neigh) {
                bond.acceptor_neigh = nb_2[iA].first;
                bond.angle_3 = angle;
            }
            if (angle < 90) {
                neighbour_distances_and_angles_are_good = false;
                break;
            }
        }
        if (!neighbour_distances_and_angles_are_good)
            break;
    }

    return std::pair<bool, h_bond_gemmi>(neighbour_distances_and_angles_are_good && good_donor_acceptor_dist, bond);
}

} // namespace h_bonds
} // namespace coot