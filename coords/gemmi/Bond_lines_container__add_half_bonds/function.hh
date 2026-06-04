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

    void add_bond_half_gemmi(int colour, const coot::Cartesian& first, const coot::Cartesian& second,
                             graphics_line_t::cylinder_class_t cc,
                             int model_number, int atom_index_1, int atom_index_2,
                             bool draw_begin, bool draw_end)
    {
        addBond(colour, first, second, cc, model_number, atom_index_1, atom_index_2, draw_begin, draw_end);
    }

    // Count total bonds across all bond groups
    int total_bond_count() const {
        int total = 0;
        for (const auto& bl : bonds) total += bl.size();
        return total;
    }
};

// gemmi port of Bond_lines_container::add_half_bonds
// Uses gemmi::CRA instead of mmdb::Atom* for atom references
void add_half_bonds_gemmi(
    Bond_lines_container_gemmi& blc,
    const coot::Cartesian &atom_1_pos,
    const coot::Cartesian &atom_2_pos,
    gemmi::CRA cra1,
    gemmi::CRA cra2,
    graphics_line_t::cylinder_class_t cc,
    int model_number,
    int atom_index_1,
    int atom_index_2,
    int atom_colour_type,
    int udd_user_defined_atom_colour_index_handle,
    my_atom_colour_map_t *atom_colour_map_p,
    bool add_begin_end_cap,
    bool add_end_end_cap);

inline void add_half_bonds_gemmi(
    Bond_lines_container_gemmi& blc,
    const coot::Cartesian &atom_1_pos,
    const coot::Cartesian &atom_2_pos,
    gemmi::CRA cra1,
    gemmi::CRA cra2,
    graphics_line_t::cylinder_class_t cc,
    int model_number,
    int atom_index_1,
    int atom_index_2,
    int atom_colour_type,
    int udd_user_defined_atom_colour_index_handle,
    my_atom_colour_map_t *atom_colour_map_p,
    bool add_begin_end_cap,
    bool add_end_end_cap) {

    coot::Cartesian bond_mid_point = atom_1_pos.mid_point(atom_2_pos);

    int col_1 = atom_colour_gemmi(cra1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
    blc.add_bond_half_gemmi(col_1, atom_1_pos, bond_mid_point, cc, model_number, atom_index_1, atom_index_2, add_begin_end_cap, false);

    int col_2 = atom_colour_gemmi(cra2, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
    blc.add_bond_half_gemmi(col_2, bond_mid_point, atom_2_pos, cc, model_number, atom_index_1, atom_index_2, false, add_end_end_cap);
}