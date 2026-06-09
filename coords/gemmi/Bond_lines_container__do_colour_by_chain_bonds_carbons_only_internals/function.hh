#pragma once

#include <gemmi/model.hpp>
#include <string>
#include <vector>

#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/Bond_lines.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__atom_colour/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__add_bond_by_dictionary_maybe/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__do_colour_by_chain_bonds_internals_goodsell_mode/gemmi/function.hh"

// Derived class to expose protected members for gemmi port
class Bond_lines_container_for_carbons : public Bond_lines_container {
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

// gemmi port of Bond_lines_container::do_colour_by_chain_bonds_carbons_only_internals
// Uses gemmi types instead of MMDB types
void do_colour_by_chain_bonds_carbons_only_internals_gemmi(
    Bond_lines_container_for_carbons& bonds,
    int imol, int imodel,
    int chain_idx,
    const gemmi::Atom *at1, const gemmi::Atom *at2,
    int iat_1, int iat_2,
    std::vector<std::pair<bool, gemmi::Residue *>> *het_residues_p,
    const std::string &element1,
    const std::string &element2,
    const coot::Cartesian &atom_1,
    const coot::Cartesian &atom_2,
    int atom_colour_type,
    int uddHnd,
    int udd_user_defined_atom_colour_index_handle);

inline void do_colour_by_chain_bonds_carbons_only_internals_gemmi(
    Bond_lines_container_for_carbons& bonds,
    int imol, int imodel,
    int chain_idx,
    const gemmi::Atom *at1, const gemmi::Atom *at2,
    int iat_1, int iat_2,
    std::vector<std::pair<bool, gemmi::Residue *>> *het_residues_p,
    const std::string &element1,
    const std::string &element2,
    const coot::Cartesian &atom_1,
    const coot::Cartesian &atom_2,
    int atom_colour_type,
    int uddHnd,
    int udd_user_defined_atom_colour_index_handle) {

    if (atom_colour_type == BC_CHAIN_GOODSELL) {
        // For goodsell, we need a separate derived class. 
        // Cast through the shared base - this works because both derive from Bond_lines_container
        Bond_lines_container_for_goodsell& goodsell_bonds = 
            static_cast<Bond_lines_container_for_goodsell&>(
                const_cast<Bond_lines_container&>(static_cast<const Bond_lines_container&>(bonds)));
        
        gemmi::Residue *res1 = nullptr;
        gemmi::Residue *res2 = nullptr;
        do_colour_by_chain_bonds_internals_goodsell_mode_gemmi(
            goodsell_bonds, imol, imodel, chain_idx,
            at1, at2, res1, res2,
            iat_1, iat_2,
            het_residues_p,
            element1, element2,
            atom_1, atom_2,
            uddHnd, udd_user_defined_atom_colour_index_handle,
            false, // have_dictionary
            false, false); // atom1_is_het, atom2_is_het
        return;
    }

    // bond_het_residue_by_dictionary = false (skip dictionary check for carbon-only bonds)
    bool bond_het_residue_by_dictionary = false;

    if (!bond_het_residue_by_dictionary) {
        if (element1 != element2) {
            bool is_H = false;
            bool draw_it = true;
            if (element1 == "H") is_H = true;
            if (element2 == "H") is_H = true;
            if (is_H) {
                double d = (atom_1 - atom_2).amplitude();
                if (d > 1.5)
                    draw_it = false;
            }

            if (draw_it) {
                coot::Cartesian bond_mid_point = atom_1.mid_point(atom_2);

                if (element1 != "C") {
                    if (element2 != "C") {
                        // half bonds, e.g. N-O, N-H, O-H
                        if (is_H) {
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            bonds.addBond(BC_COLOR_HYDROGEN_GREY, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                        } else {
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            gemmi::CRA cra1{nullptr, nullptr, const_cast<gemmi::Atom*>(at1)};
                            int non_c_col = atom_colour_gemmi(cra1, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                            bonds.bonds_size_colour_check(non_c_col);
                            bonds.addBond(non_c_col, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                            gemmi::CRA cra2{nullptr, nullptr, const_cast<gemmi::Atom*>(at2)};
                            non_c_col = atom_colour_gemmi(cra2, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                            bonds.bonds_size_colour_check(non_c_col);
                            bonds.addBond(non_c_col, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                        }
                    } else {
                        // frequent: at1 not C, at2 is C
                        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                        gemmi::CRA cra1{nullptr, nullptr, const_cast<gemmi::Atom*>(at1)};
                        int non_c_col = atom_colour_gemmi(cra1, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                        bonds.bonds_size_colour_check(non_c_col);
                        bonds.addBond(non_c_col, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                        bonds.bonds_size_colour_check(chain_idx);
                        bonds.addBond(chain_idx, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                    }
                } else {
                    // element1 is C
                    if (element2 != "C") {
                        // frequent
                        if (is_H) {
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            bonds.addBond(BC_COLOR_HYDROGEN_GREY, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                        } else {
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            bonds.bonds_size_colour_check(chain_idx);
                            bonds.addBond(chain_idx, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                            gemmi::CRA cra2{nullptr, nullptr, const_cast<gemmi::Atom*>(at2)};
                            int non_c_col = atom_colour_gemmi(cra2, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                            bonds.bonds_size_colour_check(non_c_col);
                            bonds.addBond(non_c_col, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                        }
                    } else {
                        // Both C with different elements - should not reach here
                        std::cout << "impossible" << std::endl;
                        bonds.bonds_size_colour_check(chain_idx);
                        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                        bonds.addBond(chain_idx, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                    }
                }
            }
        } else {
            // same element
            if (element1 == "C") {
                bonds.bonds_size_colour_check(chain_idx);
                graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                bonds.addBond(chain_idx, atom_1, atom_2, cc, imodel, iat_1, iat_2);
            } else {
                // same element, not carbon, not hydrogen
                if (element1 != "H") {
                    gemmi::CRA cra1{nullptr, nullptr, const_cast<gemmi::Atom*>(at1)};
                    int col = atom_colour_gemmi(cra1, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                    bonds.bonds_size_colour_check(col);
                    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                    bonds.addBond(col, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                }
            }
        }
    }
}