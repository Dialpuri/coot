#pragma once

#include <gemmi/model.hpp>
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/Bond_lines.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/Bond_lines_container__atom_colour/gemmi/function.hh"

// Helper class that inherits from Bond_lines_container to expose protected members
class Bond_lines_container_gemmi : public Bond_lines_container {
public:
    using Bond_lines_container::bonds;
    using Bond_lines_container::addBond;
};

// gemmi port of Bond_lines_container::handle_MET_or_MSE_case
// Takes CRA to carry parent context (gemmi::Atom has no parent pointer).
void handle_MET_or_MSE_case_gemmi(
    Bond_lines_container_gemmi& blc,
    gemmi::CRA cra,
    int model_number,
    int atom_colour_type,
    int udd_user_defined_atom_colour_index_handle,
    my_atom_colour_map_t* atom_colour_map_p
);

inline void handle_MET_or_MSE_case_gemmi(
    Bond_lines_container_gemmi& blc,
    gemmi::CRA cra,
    int model_number,
    int atom_colour_type,
    int udd_user_defined_atom_colour_index_handle,
    my_atom_colour_map_t* atom_colour_map_p
) {
    gemmi::Atom* atom = cra.atom;
    gemmi::Residue* residue = cra.residue;
    if (!atom || !residue) return;

    std::string atom_name(atom->name);
    std::string residue_name(residue->name);

    if (residue_name == "MET" || residue_name == "MSE" || residue_name == "MSO") {
        if (atom_name == "SE  " || atom_name == "SD" || atom_name == " SD ") {
            int col = atom_colour_gemmi(cra, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);

            for (gemmi::Atom& res_atom : residue->atoms) {
                std::string table_atom_name(res_atom.name);
                if (table_atom_name == " CG " || table_atom_name == "CG" ||
                    table_atom_name == " CE " || table_atom_name == "CE") {
                    coot::Cartesian cart_at1(atom->pos.x, atom->pos.y, atom->pos.z);
                    coot::Cartesian cart_at2(res_atom.pos.x, res_atom.pos.y, res_atom.pos.z);

                    char alt1 = atom->altloc;
                    char alt2 = res_atom.altloc;
                    if ((alt1 == 0) || (alt2 == 0) || (alt1 == alt2)) {
                        coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
                        gemmi::CRA res_cra{cra.chain, residue, &res_atom};
                        int colc = atom_colour_gemmi(res_cra, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);

                        float bond_length = (cart_at1 - cart_at2).amplitude();
                        if (bond_length < 3.0f) {
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            blc.addBond(col, cart_at1, bond_mid_point, cc, model_number, -1, -1);
                            blc.addBond(colc, bond_mid_point, cart_at2, cc, model_number, -1, -1);
                        }
                    }
                }
            }
        }
    }
    if (residue_name == "CYS") {
        int col = atom_colour_gemmi(cra, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);

        if (atom_name == " SG " || atom_name == "SG") {
            for (gemmi::Atom& res_atom : residue->atoms) {
                std::string table_atom_name(res_atom.name);
                if (table_atom_name == " CB " || table_atom_name == "CB") {
                    coot::Cartesian cart_at1(atom->pos.x, atom->pos.y, atom->pos.z);
                    coot::Cartesian cart_at2(res_atom.pos.x, res_atom.pos.y, res_atom.pos.z);

                    char alt1 = atom->altloc;
                    char alt2 = res_atom.altloc;
                    if ((alt1 == 0) || (alt2 == 0) || (alt1 == alt2)) {
                        float len2 = (cart_at1 - cart_at2).amplitude_squared();
                        if (len2 < 16.0f) {
                            coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
                            gemmi::CRA res_cra{cra.chain, residue, &res_atom};
                            int colc = atom_colour_gemmi(res_cra, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
                            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                            blc.addBond(col, cart_at1, bond_mid_point, cc, model_number, -1, -1);
                            blc.addBond(colc, bond_mid_point, cart_at2, cc, model_number, -1, -1);
                        }
                    }
                }
            }
        }
    }
}
