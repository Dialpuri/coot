#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "clipper/core/coords.h"
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/Bond_lines.hh"
#include "geometry/protein-geometry.hh"

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__get_neighb_normal/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__atom_colour/gemmi/function.hh"

/// Helper: add a dashed bond to bonds_storage
/// Mirrors Bond_lines_container::add_dashed_bond but operates on the
/// bonds_storage map directly (no Bond_lines_container needed).
static void add_dashed_bond_to_storage(
    std::map<int, Bond_lines>& bonds_storage,
    int col,
    const coot::Cartesian& cart_start,
    const coot::Cartesian& cart_end,
    int half_bond_type,
    graphics_line_t::cylinder_class_t cc,
    int imodel,
    int idx_1_mol_indexing,
    int idx_2_mol_indexing)
{
    clipper::Coord_orth pos_1(cart_start.x(), cart_start.y(), cart_start.z());
    clipper::Coord_orth pos_2(cart_end.x(), cart_end.y(), cart_end.z());
    clipper::Coord_orth v(pos_2.x() - pos_1.x(), pos_2.y() - pos_1.y(), pos_2.z() - pos_1.z());
    double dist = std::sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z());
    int n_segs = static_cast<int>(dist / 0.5) + 2;
    if (n_segs < 2) n_segs = 2;
    double inv_n = 1.0 / n_segs;
    clipper::Coord_orth step(v.x() * inv_n, v.y() * inv_n, v.z() * inv_n);
    clipper::Coord_orth current = pos_1;
    for (int seg = 0; seg < n_segs; seg++) {
        clipper::Coord_orth next_pos(current.x() + step.x(), current.y() + step.y(), current.z() + step.z());
        bonds_storage[col].add_bond(
            coot::CartesianPair(coot::Cartesian(current), coot::Cartesian(next_pos)),
            cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
        current = next_pos;
    }
}

/// gemmi port of Bond_lines_container::add_double_bond
/// Adds 2 (same-element) or 4 (different-element) bond graphics_line_t objects
/// to represent a double bond between two atoms.
///
/// Parameters:
///   bonds_storage — map from colour index to Bond_lines vector (mirrors Bond_lines_container::bonds)
///   imol, imodel  — molecule and model indices
///   iat_1, iat_2  — atom indices within the atoms vector
///   atoms         — gemmi atoms for the residue
///   res_name      — residue name for geometry lookup
///   atom_colour_type — colouring scheme
///   atom_colour_map_p — user-defined colour map (nullable)
///   udd_user_defined_atom_colour_index_handle — UDD colour index
///   geom          — protein geometry for neighbour normal computation
///   is_deloc      — whether the bond involves delocalized electrons
inline void add_double_bond_gemmi(
    std::map<int, Bond_lines>& bonds_storage,
    int imol,
    int imodel,
    int iat_1,
    int iat_2,
    const std::vector<gemmi::Atom>& atoms,
    const std::string& res_name,
    int atom_colour_type,
    my_atom_colour_map_t *atom_colour_map_p,
    int udd_user_defined_atom_colour_index_handle,
    coot::protein_geometry* geom,
    bool is_deloc) {

    std::string ele_1 = atoms[iat_1].element.name();
    std::string ele_2 = atoms[iat_2].element.name();

    graphics_line_t::cylinder_class_t cc = graphics_line_t::DOUBLE;

    int idx_1_mol_indexing = iat_1;
    int idx_2_mol_indexing = iat_2;

    try {

        // perp_n is the direction of the offset (from the atom position) of the start and
        // finish points in the plane of the double bond.

        clipper::Coord_orth pos_at_1(atoms[iat_1].pos.x, atoms[iat_1].pos.y, atoms[iat_1].pos.z);
        clipper::Coord_orth pos_at_2(atoms[iat_2].pos.x, atoms[iat_2].pos.y, atoms[iat_2].pos.z);

        // Compute neighbour normal
        clipper::Coord_orth n_n(0, 0, 1); // default fallback
        bool also_2nd_order = true;
        if (geom != nullptr) {
            n_n = get_neighb_normal_gemmi(*geom, imol, iat_1, iat_2, atoms, res_name, also_2nd_order);
        }

        clipper::Coord_orth b(pos_at_1.x() - pos_at_2.x(), pos_at_1.y() - pos_at_2.y(), pos_at_1.z() - pos_at_2.z());
        double b_len = std::sqrt(b.x()*b.x() + b.y()*b.y() + b.z()*b.z());
        double inv_b = 1.0 / b_len;
        clipper::Coord_orth b_n(b.x() * inv_b, b.y() * inv_b, b.z() * inv_b);
        clipper::Coord_orth perp_n(clipper::Coord_orth::cross(n_n, b_n));

        // Build CRA for atom_colour_gemmi call
        gemmi::CRA cra_1;
        cra_1.atom = const_cast<gemmi::Atom*>(&atoms[iat_1]);
        cra_1.residue = nullptr;
        cra_1.chain = nullptr;

        int col = atom_colour_gemmi(cra_1, atom_colour_type,
                                    udd_user_defined_atom_colour_index_handle,
                                    atom_colour_map_p);

        double offset = 0.13;
        clipper::Coord_orth pt_1_1(pos_at_1.x() - offset * perp_n.x(),
                                   pos_at_1.y() - offset * perp_n.y(),
                                   pos_at_1.z() - offset * perp_n.z());
        clipper::Coord_orth pt_1_2(pos_at_1.x() + offset * perp_n.x(),
                                   pos_at_1.y() + offset * perp_n.y(),
                                   pos_at_1.z() + offset * perp_n.z());
        clipper::Coord_orth pt_2_1(pos_at_2.x() - offset * perp_n.x(),
                                   pos_at_2.y() - offset * perp_n.y(),
                                   pos_at_2.z() - offset * perp_n.z());
        clipper::Coord_orth pt_2_2(pos_at_2.x() + offset * perp_n.x(),
                                   pos_at_2.y() + offset * perp_n.y(),
                                   pos_at_2.z() + offset * perp_n.z());

        if (ele_1 == ele_2) {
            // simple double bond (e.g. C=C)
            bonds_storage[col].add_bond(
                coot::CartesianPair(coot::Cartesian(pt_1_1), coot::Cartesian(pt_2_1)),
                cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            if (!is_deloc)
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2)),
                    cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            else
                add_dashed_bond_to_storage(bonds_storage, col,
                    coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2),
                    0, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
        } else {
            // different elements — draw double half bonds, e.g. C=O
            clipper::Coord_orth bond_mid_point(
                0.5 * (pos_at_1.x() + pos_at_2.x()),
                0.5 * (pos_at_1.y() + pos_at_2.y()),
                0.5 * (pos_at_1.z() + pos_at_2.z()));
            clipper::Coord_orth mp_1(bond_mid_point.x() - offset * perp_n.x(),
                                     bond_mid_point.y() - offset * perp_n.y(),
                                     bond_mid_point.z() - offset * perp_n.z());
            clipper::Coord_orth mp_2(bond_mid_point.x() + offset * perp_n.x(),
                                     bond_mid_point.y() + offset * perp_n.y(),
                                     bond_mid_point.z() + offset * perp_n.z());
            if (!is_deloc) {
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_1_1), coot::Cartesian(mp_1)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_1_2), coot::Cartesian(mp_2)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);

                gemmi::CRA cra_2;
                cra_2.atom = const_cast<gemmi::Atom*>(&atoms[iat_2]);
                cra_2.residue = nullptr;
                cra_2.chain = nullptr;

                col = atom_colour_gemmi(cra_2, atom_colour_type,
                                        udd_user_defined_atom_colour_index_handle,
                                        atom_colour_map_p);
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_2_1), coot::Cartesian(mp_1)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_2_2), coot::Cartesian(mp_2)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            } else {
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_1_1), coot::Cartesian(mp_1)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);

                add_dashed_bond_to_storage(bonds_storage, col,
                    coot::Cartesian(pt_1_2), coot::Cartesian(mp_2),
                    1, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);

                gemmi::CRA cra_2;
                cra_2.atom = const_cast<gemmi::Atom*>(&atoms[iat_2]);
                cra_2.residue = nullptr;
                cra_2.chain = nullptr;

                col = atom_colour_gemmi(cra_2, atom_colour_type,
                                        udd_user_defined_atom_colour_index_handle,
                                        atom_colour_map_p);
                bonds_storage[col].add_bond(
                    coot::CartesianPair(coot::Cartesian(pt_2_1), coot::Cartesian(mp_1)),
                    cc, true, false, imodel, idx_1_mol_indexing, idx_2_mol_indexing);

                add_dashed_bond_to_storage(bonds_storage, col,
                    coot::Cartesian(pt_2_2), coot::Cartesian(mp_2),
                    2, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            }
        }

    } catch (const std::runtime_error &rte) {
        std::cout << "caught exception add_double_bond_gemmi(): " << rte.what() << std::endl;
    }
}
