#pragma once

#include "mini-mol/mini-mol.hh"
#include <gemmi/math.hpp>
#include <vector>
#include <string>
#include <utility>

namespace coot {
namespace minimol {

// gemmi port of coot::minimol::molecule::molecule constructor
// Input uses gemmi::Vec3 instead of clipper::Coord_orth
inline molecule molecule_gemmi(
    const std::vector<std::pair<gemmi::Vec3, float>> &atom_list_with_estimated_b_factors,
    const std::string& residue_type,
    std::string atom_name,
    std::string chain_id,
    const std::string &ele)
{
    molecule mol;

    // Constructing a fragment from a chain_id sets residues_offset to 0
    // but doesn't add any residues
    mol.fragments.push_back(fragment(chain_id));
    std::string element = ele;

    // Each atom goes in its own residue (residue number offset by one
    // c.f. the atom vector index)
    for (unsigned int i = 0; i < atom_list_with_estimated_b_factors.size(); i++) {
        auto b_est = atom_list_with_estimated_b_factors[i].second;
        mol.fragments[0][i+1] = residue(i+1); // atoms start at 0, residues at 1.
        mol.fragments[0][i+1].name = residue_type;  // not "WAT" says EJD - 030624
        mol.fragments[0][i+1].addatom(
            atom_name,
            element,
            atom_list_with_estimated_b_factors[i].first.x,
            atom_list_with_estimated_b_factors[i].first.y,
            atom_list_with_estimated_b_factors[i].first.z,
            std::string(""),
            1.0f,
            b_est); // pass this? 20090201
    }
    // have_cell and have_spacegroup default to 0; skip explicit private-member writes

    return mol;
}

} // namespace minimol
} // namespace coot
