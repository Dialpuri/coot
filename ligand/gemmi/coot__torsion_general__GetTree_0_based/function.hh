#pragma once

#include <vector>
#include <cmath>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "ccp4mg-utils/mgtree.h"

namespace coot {
namespace torsion_general {

static std::vector<std::vector<int>>
compute_contact_indices_gemmi(const gemmi::Residue& res) {
    int n = static_cast<int>(res.atoms.size());
    std::vector<std::vector<int>> v(n);

    // Distance thresholds from original get_contact_indices()
    double max_dist = 1.9;
    double max_dist_bond_to_H = 1.42;

    // Check for MSE (selenomethionine)
    if (res.name == "MSE")
        max_dist = 2.0;

    // Separate H and non-H atoms (keeping original indices)
    std::vector<int> non_H_indices;
    std::vector<int> H_indices;
    for (int i = 0; i < n; i++) {
        if (res.atoms[i].element.is_hydrogen()) {
            H_indices.push_back(i);
        } else {
            non_H_indices.push_back(i);
        }
    }

    // Non-H to non-H contacts (within same residue)
    for (int ii = 0; ii < static_cast<int>(non_H_indices.size()); ii++) {
        int i = non_H_indices[ii];
        for (int jj = ii + 1; jj < static_cast<int>(non_H_indices.size()); jj++) {
            int j = non_H_indices[jj];
            double dx = res.atoms[i].pos.x - res.atoms[j].pos.x;
            double dy = res.atoms[i].pos.y - res.atoms[j].pos.y;
            double dz = res.atoms[i].pos.z - res.atoms[j].pos.z;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist >= 0.1 && dist <= max_dist) {
                v[i].push_back(j);
            }
        }
    }

    // H to non-H contacts
    for (int ii = 0; ii < static_cast<int>(H_indices.size()); ii++) {
        int i = H_indices[ii];
        for (int jj = 0; jj < static_cast<int>(non_H_indices.size()); jj++) {
            int j = non_H_indices[jj];
            double dx = res.atoms[i].pos.x - res.atoms[j].pos.x;
            double dy = res.atoms[i].pos.y - res.atoms[j].pos.y;
            double dz = res.atoms[i].pos.z - res.atoms[j].pos.z;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist >= 0.1 && dist <= max_dist_bond_to_H) {
                v[i].push_back(j);
            }
        }
    }

    return v;
}

// Gemmi port of coot::torsion_general::GetTree_0_based()
// Extracts all residue atom coordinates into a Tree with 0-based indexing.
inline Tree GetTree_0_based_gemmi(gemmi::Residue* residue_p)
{
    Tree tree;
    if (residue_p) {
        std::vector< ::Cartesian > coords;
        for (const auto& atom : residue_p->atoms) {
            coords.emplace_back(atom.pos.x, atom.pos.y, atom.pos.z);
        }
        auto contact_indices = compute_contact_indices_gemmi(*residue_p);
        int base_index = 0;
        tree.SetCoords(coords, base_index, contact_indices);
    }
    return tree;
}

} // namespace torsion_general
} // namespace coot
