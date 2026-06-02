#pragma once

#include <gemmi/math.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

struct TreeVertex {
    gemmi::Vec3 coord;
    TreeVertex* parent;
    std::vector<TreeVertex*> children;
    int id;
    double parent_dist;
    double parent_bond_angle;
    double parent_dihedral_angle;

    TreeVertex() : parent(nullptr), id(0), parent_dist(0.0),
                   parent_bond_angle(0.0), parent_dihedral_angle(0.0) {}
};

struct Tree {
    std::vector<TreeVertex*> coords;
};

// Forward declaration
inline std::ostream& print_vertex_gemmi(std::ostream& c, const TreeVertex& a);

inline std::ostream& print_tree_gemmi(std::ostream& c, const Tree& a) {
    c << std::setprecision(6);
    c << std::fixed;

    for (unsigned int i = 0; i < a.coords.size(); i++) {
        const gemmi::Vec3& v = a.coords[i]->coord;
        if (a.coords[i]->parent != nullptr) {
            const gemmi::Vec3& pv = a.coords[i]->parent->coord;
            c << "Atom[" << i << "]: "
              << v.x << " " << v.y << " " << v.z
              << " has parent atom[" << a.coords[i]->parent->id << "]: "
              << pv.x << " " << pv.y << " " << pv.z << "\n";
        } else {
            c << "Atom[" << i << "]: "
              << v.x << " " << v.y << " " << v.z
              << " has no parent\n";
        }
    }

    c << "\nTrees ....\n";
    for (unsigned int i = 0; i < a.coords.size(); i++) {
        if (a.coords[i]->parent == nullptr) {
            print_vertex_gemmi(c, *a.coords[i]);
        }
    }

    c.unsetf(std::ios::fixed | std::ios::scientific);
    return c;
}

inline std::ostream& print_vertex_gemmi(std::ostream& c, const TreeVertex& a) {
    c << std::setprecision(6);
    c << std::fixed;
    constexpr double rtod = 180.0 / M_PI;

    for (unsigned int i = 0; i < a.children.size(); i++) {
        print_vertex_gemmi(c, *a.children[i]);
    }
    c << "\n";

    int depth = 0;
    const TreeVertex* p = a.parent;
    while (p != nullptr) {
        ++depth;
        p = p->parent;
    }
    for (int i = 0; i < depth; ++i) {
        c << " ";
    }
    c << std::setw(3) << a.id << ": [Parent "
      << (a.parent ? a.parent->id : -1)
      << "] Length:" << a.parent_dist
      << " Angle:" << a.parent_bond_angle * rtod
      << " Torsion:" << a.parent_dihedral_angle * rtod << "\n";

    for (unsigned int i = 0; i < a.children.size(); i++) {
        print_vertex_gemmi(c, *a.children[i]);
    }

    c.unsetf(std::ios::fixed | std::ios::scientific);
    return c;
}