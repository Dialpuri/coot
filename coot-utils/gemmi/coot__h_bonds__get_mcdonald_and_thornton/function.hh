#pragma once

#include <gemmi/model.hpp>
#include <vector>

#include "geometry/protein-geometry.hh"

namespace coot {

// Unified h_bond struct for gemmi port (has all fields the test needs)
struct mcdonald_h_bond {
    double dist;
    double angle_1;
    double angle_2;
    double angle_3;

    const gemmi::Atom *hb_hydrogen;
    const gemmi::Atom *acceptor;
    const gemmi::Atom *donor;
    const gemmi::Atom *donor_neigh;
    const gemmi::Atom *acceptor_neigh;

    bool ligand_atom_is_donor;
    bool bond_has_hydrogen_flag;

    mcdonald_h_bond()
        : dist(0), angle_1(0), angle_2(-1), angle_3(-1),
          hb_hydrogen(nullptr), acceptor(nullptr), donor(nullptr),
          donor_neigh(nullptr), acceptor_neigh(nullptr),
          ligand_atom_is_donor(false), bond_has_hydrogen_flag(false) {}
};

std::vector<mcdonald_h_bond>
get_mcdonald_and_thornton_gemmi(
    const std::vector<gemmi::CRA>& sel_1_atoms,
    const std::vector<gemmi::CRA>& sel_2_atoms,
    gemmi::Model& model,
    const coot::protein_geometry& geom,
    int imol,
    double max_dist);

} // namespace coot
