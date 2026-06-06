#pragma once

#include <gemmi/model.hpp>
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/Bond_lines.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__add_bond_by_dictionary_maybe/gemmi/function.hh"
#include "geometry/protein-geometry.hh"
#include <vector>
#include <string>
#include <utility>

// Derived class to expose protected members for gemmi port
class Bond_lines_container_for_goodsell : public Bond_lines_container {
public:
    using Bond_lines_container::bonds;
    using Bond_lines_container::bonds_size_colour_check;
    using Bond_lines_container::addBond;

    int total_bond_count() const {
        int total = 0;
        for (const auto& bl : bonds) total += bl.size();
        return total;
    }
};

// gemmi port of Bond_lines_container::do_colour_by_chain_bonds_internals_goodsell_mode
// Uses gemmi types instead of MMDB types
// Note: have_dictionary is passed as parameter (it's private in Bond_lines_container)
void do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
    Bond_lines_container_for_goodsell& bonds,
    int imol, int imodel,
    int chain_idx,
    const gemmi::Atom* atom1, const gemmi::Atom* atom2,
    gemmi::Residue* residue_1, gemmi::Residue* residue_2,
    int iat_1, int iat_2,
    std::vector<std::pair<bool, gemmi::Residue*>>* het_residues_p,
    const std::string& element1,
    const std::string& element2,
    const coot::Cartesian& atom_pos_1,
    const coot::Cartesian& atom_pos_2,
    int uddHnd,
    int udd_user_defined_atom_colour_index_handle,
    bool have_dictionary,
    bool atom1_is_het = false,
    bool atom2_is_het = false);

inline void do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
    Bond_lines_container_for_goodsell& bonds,
    int imol, int imodel,
    int chain_idx,
    const gemmi::Atom* atom1, const gemmi::Atom* atom2,
    gemmi::Residue* residue_1, gemmi::Residue* residue_2,
    int iat_1, int iat_2,
    std::vector<std::pair<bool, gemmi::Residue*>>* het_residues_p,
    const std::string& element1,
    const std::string& element2,
    const coot::Cartesian& atom_pos_1,
    const coot::Cartesian& atom_pos_2,
    int /*uddHnd*/,
    int /*udd_user_defined_atom_colour_index_handle*/,
    bool have_dictionary,
    bool atom1_is_het,
    bool atom2_is_het) {

    coot::protein_geometry geom;
    geom.init_standard();

    bool bond_het_residue_by_dictionary =
        add_bond_by_dictionary_maybe_gemmi(
            imol, atom1, atom2,
            residue_1, residue_2, het_residues_p,
            &geom, have_dictionary,
            atom1_is_het, atom2_is_het);

    if (!bond_het_residue_by_dictionary) {
        bool draw_it = true;
        bool is_H = false;
        if (element1 == " H") is_H = true;
        if (element2 == " H") is_H = true;
        if (is_H) {
            double d = (atom_pos_1 - atom_pos_2).amplitude();
            if (d > 1.5)
                draw_it = false;
        }

        if (draw_it) {
            coot::Cartesian bond_mid_point = atom_pos_1.mid_point(atom_pos_2);
            if (element1 != " C") {
                if (element2 != " C") {
                    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                    int non_c_col = 2 * chain_idx + 1;
                    bonds.bonds_size_colour_check(non_c_col);
                    bonds.addBond(non_c_col, atom_pos_1, atom_pos_2, cc, imodel, iat_1, iat_2);
                } else {
                    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                    int non_c_col = 2 * chain_idx + 1;
                    bonds.bonds_size_colour_check(non_c_col);
                    bonds.addBond(non_c_col, atom_pos_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                    int c_col = 2 * chain_idx;
                    bonds.bonds_size_colour_check(c_col);
                    bonds.addBond(c_col, atom_pos_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                }
            } else {
                // first atom *was* carbon
                if (element2 == " C") {
                    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                    int c_col = 2 * chain_idx;
                    bonds.bonds_size_colour_check(c_col);
                    bonds.addBond(c_col, atom_pos_1, atom_pos_2, cc, imodel, iat_1, iat_2);
                } else {
                    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                    int non_c_col = 2 * chain_idx + 1;
                    bonds.bonds_size_colour_check(non_c_col);
                    bonds.addBond(non_c_col, atom_pos_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                    int c_col = 2 * chain_idx;
                    bonds.bonds_size_colour_check(c_col);
                    bonds.addBond(c_col, atom_pos_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                }
            }
            // PutUDData: NO_EQUIVALENT in gemmi; skipped.
        }
    }
}