#pragma once

#include <string>
#include <vector>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__distance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__angle/gemmi/function.hh"

namespace coot {
namespace h_bonds_gemmi {

class h_bond_gemmi {
public:
    const gemmi::Atom* hb_hydrogen; // McDonald and Thornton H-bond algorithm
    const gemmi::Atom* donor;
    const gemmi::Atom* acceptor;
    const gemmi::Atom* donor_neigh;
    const gemmi::Atom* acceptor_neigh;
    double angle_1;  // degrees
    double angle_2;
    double angle_3;
    double dist;  // H-bond length
    bool ligand_atom_is_donor; // for use when hb_hydrogen is NULL -
                               // no hydrogens in H-bond analysis.
    bool hydrogen_is_ligand_atom;
    bool bond_has_hydrogen_flag;
    
    h_bond_gemmi() {
        hb_hydrogen = nullptr;
        donor = nullptr;
        acceptor = nullptr;
        donor_neigh = nullptr;
        acceptor_neigh = nullptr;
        ligand_atom_is_donor = 1;
        angle_1 = -1;
        angle_2 = -1;
        angle_3 = -1;
        hydrogen_is_ligand_atom = 0; // no hydrogen
        bond_has_hydrogen_flag = 0;
    }
    
    h_bond_gemmi(const gemmi::Atom* d, const gemmi::Atom* a) {
        hb_hydrogen = nullptr;
        donor = d;
        acceptor = a;
        donor_neigh = nullptr;
        acceptor_neigh = nullptr;
        ligand_atom_is_donor = 0;
        dist = -1;
        angle_1 = -1;
        angle_2 = -1;
        angle_3 = -1; 
        hydrogen_is_ligand_atom = 0; // no hydrogen
        bond_has_hydrogen_flag = 0;
    }

    // for McDonald and Thornton H-bonds, where the Hs are explicit.
    // One first of these atoms is the hydrogen, the other is the
    // acceptor.
    //
    // pass ligand_atom_is_H_flag as 1 when ligand atom is the H.
    // 
    h_bond_gemmi(const gemmi::Atom* h, const gemmi::Atom* a, bool ligand_atom_is_H_flag) {
        hb_hydrogen = h;
        bond_has_hydrogen_flag = 1;
        acceptor = a;
        donor = nullptr;
        donor_neigh = nullptr;
        acceptor_neigh = nullptr;
        hydrogen_is_ligand_atom = ligand_atom_is_H_flag;
        ligand_atom_is_donor = ligand_atom_is_H_flag;
        dist = -1;
        angle_1 = -1; 
        angle_2 = -1; 
        angle_3 = -1; 
    }
};

// Ported version of make_h_bond_from_environment_residue_hydrogen
// Takes gemmi atoms instead of mmdb atoms
std::pair<bool, h_bond_gemmi> 
make_h_bond_from_environment_residue_hydrogen_gemmi(const gemmi::Atom* at_1, // acceptor on ligand
                                                    const gemmi::Atom* at_2, // H on residue
                                                    const std::vector<std::pair<const gemmi::Atom*, float>>& nb_1,
                                                    const std::vector<std::pair<const gemmi::Atom*, float>>& nb_2);

} // namespace h_bonds_gemmi
} // namespace coot