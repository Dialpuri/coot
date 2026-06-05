#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include "ccp4mg-utils/mgtree.h"

namespace coot { namespace torsion_general {

inline std::vector<std::vector<int>>
get_contact_indices_gemmi(const gemmi::Residue& residue) {
    std::vector<std::vector<int>> v;
    float min_dist = 0.1;
    float max_dist = 1.9;
    float max_dist_bond_to_H = 1.42;
    if (residue.name == "MSE")
        max_dist = 2.0;

    std::vector<std::pair<size_t, gemmi::Position>> h_atoms;
    std::vector<std::pair<size_t, gemmi::Position>> nh_atoms;
    
    for (size_t i = 0; i < residue.atoms.size(); ++i) {
        const auto& atom = residue.atoms[i];
        if (atom.element.is_hydrogen()) {
            h_atoms.emplace_back(i, atom.pos);
        } else {
            nh_atoms.emplace_back(i, atom.pos);
        }
    }

    v.resize(residue.atoms.size());

    for (size_t i = 0; i < nh_atoms.size(); ++i) {
        for (size_t j = i + 1; j < nh_atoms.size(); ++j) {
            double dist_sq = nh_atoms[i].second.dist_sq(nh_atoms[j].second);
            if (dist_sq >= min_dist * min_dist && dist_sq <= max_dist * max_dist) {
                v[nh_atoms[i].first].push_back(nh_atoms[j].first);
                v[nh_atoms[j].first].push_back(nh_atoms[i].first);
            }
        }
    }

    for (const auto& h : h_atoms) {
        for (const auto& nh : nh_atoms) {
            double dist_sq = h.second.dist_sq(nh.second);
            if (dist_sq >= min_dist * min_dist && dist_sq <= max_dist_bond_to_H * max_dist_bond_to_H) {
                v[h.first].push_back(nh.first);
            }
        }
    }

    return v;
}

inline Tree GetTree_gemmi(const gemmi::Residue& residue) {
    Tree tree;
    std::vector<::Cartesian> coords;
    for (const auto& atom : residue.atoms) {
        coords.emplace_back(atom.pos.x, atom.pos.y, atom.pos.z);
    }
    auto contact_indices = get_contact_indices_gemmi(residue);
    tree.SetCoords(coords, 0, contact_indices);
    return tree;
}

}}