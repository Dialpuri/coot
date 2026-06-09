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

/// gemmi port of Bond_lines_container::add_triple_bond
/// Adds 6 bond graphics_line_t objects to represent a triple bond between two atoms.
/// When elements match: 3 parallel lines from atom1 to atom2.
/// When elements differ: 3 lines from atom1 to midpoint, 3 lines from atom2 to midpoint.
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
void add_triple_bond_gemmi(
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
    coot::protein_geometry* geom);

inline void add_triple_bond_gemmi(
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
    coot::protein_geometry* geom) {

    graphics_line_t::cylinder_class_t cc = graphics_line_t::TRIPLE;

    // Get element names using gemmi (proper string comparison, unlike MMDB's pointer comparison)
    std::string ele_1 = atoms[iat_1].element.name();
    std::string ele_2 = atoms[iat_2].element.name();

    try {

        int idx_1_mol_indexing = iat_1;
        int idx_2_mol_indexing = iat_2;

        bool also_2nd_order = true; // linear nature of bonds in triple bond

        clipper::Coord_orth pos_at_1(atoms[iat_1].pos.x, atoms[iat_1].pos.y, atoms[iat_1].pos.z);
        clipper::Coord_orth pos_at_2(atoms[iat_2].pos.x, atoms[iat_2].pos.y, atoms[iat_2].pos.z);

        // Compute neighbour normal using the existing gemmi port
        clipper::Coord_orth n_n = clipper::Coord_orth(0, 0, 1); // default fallback
        if (geom != nullptr) {
            n_n = get_neighb_normal_gemmi(*geom, imol, iat_1, iat_2, atoms, res_name, also_2nd_order);
        }

        clipper::Coord_orth b(pos_at_1 - pos_at_2);
        clipper::Coord_orth b_n(b.unit());
        clipper::Coord_orth perp_n(clipper::Coord_orth::cross(n_n, b_n));

        // Build CRA for atom_colour_gemmi call
        gemmi::CRA cra_1;
        cra_1.atom = const_cast<gemmi::Atom*>(&atoms[iat_1]);
        cra_1.residue = nullptr;
        cra_1.chain = nullptr;

        int col = atom_colour_gemmi(cra_1, atom_colour_type,
                                    udd_user_defined_atom_colour_index_handle,
                                    atom_colour_map_p);
        double offset = 0.08;

        clipper::Coord_orth pt_1_1 = pos_at_1 - offset * perp_n;
        clipper::Coord_orth pt_1_2 = pos_at_1;
        clipper::Coord_orth pt_1_3 = pos_at_1 + offset * perp_n;
        clipper::Coord_orth pt_2_1 = pos_at_2 - offset * perp_n;
        clipper::Coord_orth pt_2_2 = pos_at_2;
        clipper::Coord_orth pt_2_3 = pos_at_2 + offset * perp_n;

        if (ele_1 == ele_2) {
            // Same element (e.g. -C#C-)
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_1), coot::Cartesian(pt_2_1)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_3), coot::Cartesian(pt_2_3)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
        } else {
            // Different elements (e.g. -C#N)
            clipper::Coord_orth bond_mid_point = 0.5 * clipper::Coord_orth(pos_at_1 + pos_at_2);
            clipper::Coord_orth mp_1 = bond_mid_point - offset * perp_n;
            clipper::Coord_orth mp_2 = bond_mid_point;
            clipper::Coord_orth mp_3 = bond_mid_point + offset * perp_n;

            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_1), coot::Cartesian(mp_1)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_2), coot::Cartesian(mp_2)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_1_3), coot::Cartesian(mp_3)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);

            // Second colour for atom 2
            gemmi::CRA cra_2;
            cra_2.atom = const_cast<gemmi::Atom*>(&atoms[iat_2]);
            cra_2.residue = nullptr;
            cra_2.chain = nullptr;

            col = atom_colour_gemmi(cra_2, atom_colour_type,
                                    udd_user_defined_atom_colour_index_handle,
                                    atom_colour_map_p);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_2_1), coot::Cartesian(mp_1)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_2_2), coot::Cartesian(mp_2)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            bonds_storage[col].add_bond(coot::CartesianPair(coot::Cartesian(pt_2_3), coot::Cartesian(mp_3)), cc, true, true, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
        }

    } catch (const std::runtime_error &rte) {
        std::cout << "caught exception add_triple_bond_gemmi(): " << rte.what() << std::endl;
    }
}
