#pragma once

#include <vector>
#include <string>
#include <gemmi/math.hpp>
#include "mini-mol/mini-mol.hh"

namespace coot::minimol {

inline molecule molecule_gemmi(const std::vector<gemmi::Vec3> &atom_list,
                               const std::string& residue_type,
                               std::string atom_name,
                               std::string chain_id,
                               const std::string &ele) {
    molecule mol;

    // Constructing a fragment from a chain_id sets residues_offset to 0
    // but doesn't add any residues
    mol.fragments.push_back(fragment(chain_id));
    std::string element = ele;

    // Each atom goes in its own residue (residue number offset by one
    // c.f. the atom vector index)
    for (unsigned int i = 0; i < atom_list.size(); i++) {
        mol.fragments[0][i + 1] = residue(i + 1); // atoms start at 0, residues at 1.
        mol.fragments[0][i + 1].name = residue_type; // not "WAT" says EJD - 030624
        mol.fragments[0][i + 1].addatom(atom_name,
                                        element,
                                        atom_list[i].x,
                                        atom_list[i].y,
                                        atom_list[i].z,
                                        std::string(""),
                                        1.0,
                                        30.0); // pass this? 20090201
    }
    return mol;
}

} // namespace coot::minimol
