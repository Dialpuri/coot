#pragma once

#include <gemmi/model.hpp>
#include "coords/Cartesian.hh"
#include "coords/graphics-line.hh"
#include "coords/Bond_lines.hh"

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

    int total_bond_count() const {
        int total = 0;
        for (const auto& bl : bonds) total += bl.size();
        return total;
    }

    /**
     * gemmi port of Bond_lines_container::handle_long_bonded_atom
     *
     * For atoms that don't have enough standard bonds, this searches
     * within the same residue for close atoms to bond with.
     *
     * Uses gemmi::CRA instead of mmdb::PAtom for parent context.
     *
     * @param cra        The atom and its parent chain/residue context
     * @param model_number  Model number (0-based)
     * @param udd_handle_bond              Unused in gemmi (MMDB UDD marker)
     * @param udd_handle_atom_index        Unused in gemmi (MMDB UDD index)
     * @param udd_user_defined_atom_colour_index_handle  Passed to atom_colour
     * @param atom_colour_type             Colouring scheme (COLOUR_BY_CHAIN=0, etc.)
     * @param atom_colour_map_p            Colour map
     * @return Number of bonds added (0, 2 per neighbor, or 3 for star)
     */
    int handle_long_bonded_atom_gemmi(
        gemmi::CRA cra,
        int model_number,
        int udd_handle_bond,
        int udd_handle_atom_index,
        int udd_user_defined_atom_colour_index_handle,
        int atom_colour_type,
        coot::my_atom_colour_map_t *atom_colour_map_p);
};

// Inline implementation
inline int Bond_lines_container_gemmi::handle_long_bonded_atom_gemmi(
    gemmi::CRA cra,
    int model_number,
    int /*udd_handle_bond*/,
    int /*udd_handle_atom_index*/,
    int /*udd_user_defined_atom_colour_index_handle*/,
    int atom_colour_type,
    coot::my_atom_colour_map_t *atom_colour_map_p)
{
    if (!cra.atom || !cra.residue || !cra.chain)
        return 0;

    float bond_limit = 2.16f;

    std::string atom_name(cra.atom->name);
    std::string element(cra.atom->element.name());

    // Adjust bond limit for special elements
    if (atom_name == "AS  ")
        bond_limit = 2.4f;
    if (element == "AU")
        bond_limit = 2.4f;
    if (element == "AS")
        bond_limit = 2.4f;
    if (element == "HG")
        bond_limit = 2.4f;
    if (element == "MO")
        bond_limit = 2.55f;
    if (element == "Mo")
        bond_limit = 2.55f;
    if (element == " I")
        bond_limit = 2.3f;

    float bl2 = bond_limit * bond_limit;
    float h_bl2 = 1.8f * 1.8f;

    coot::Cartesian atom_pos(cra.atom->pos.x, cra.atom->pos.y, cra.atom->pos.z);
    int bonds_added = 0;

    for (gemmi::Atom& res_atom : cra.residue->atoms) {
        if (&res_atom == cra.atom)
            continue;

        coot::Cartesian res_atom_pos(res_atom.pos.x, res_atom.pos.y, res_atom.pos.z);
        float len2 = (atom_pos - res_atom_pos).amplitude_squared();

        bool res_atom_is_hydrogen = res_atom.element.is_hydrogen();
        bool distance_ok = res_atom_is_hydrogen
                            ? (len2 < h_bl2)
                            : (len2 < bl2);

        if (!distance_ok)
            continue;

        // Check altloc compatibility: gemmi uses '\0' for no altloc (like MMDB uses "")
        char alt1 = cra.atom->altloc;
        char alt2 = res_atom.altloc;
        bool alt_ok = (alt1 == '\0' || alt2 == '\0' || alt1 == alt2);
        if (!alt_ok)
            continue;

        coot::Cartesian bond_mid_point = atom_pos.mid_point(res_atom_pos);
        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

        gemmi::CRA cra_res_atom{cra.chain, cra.residue, &res_atom};

        // Atom colour lookup — COLOUR_BY_CHAIN = 0
        // This mirrors Bond_lines_container::atom_colour for atom_colour_type=0
        int col = 0;
        if (atom_colour_type == 0 && atom_colour_map_p && cra.chain) {
            col = atom_colour_map_p->index_for_chain(cra.chain->name);
        }
        int colc = 0;
        if (atom_colour_type == 0 && atom_colour_map_p && cra_res_atom.chain) {
            colc = atom_colour_map_p->index_for_chain(cra_res_atom.chain->name);
        }

        // Add two half-bonds: atom→midpoint and midpoint→res_atom
        add_bond_half_gemmi(col, atom_pos, bond_mid_point, cc, model_number, -1, -1, false, false);
        add_bond_half_gemmi(colc, bond_mid_point, res_atom_pos, cc, model_number, -1, -1, false, false);

        bonds_added += 2;
    }

    if (bonds_added == 0) {
        // Draw "star" for unbonded atom (3 crossing lines)
        float star_size = 0.3f;
        coot::Cartesian vec_x(star_size, 0.0f, 0.0f);
        coot::Cartesian vec_y(0.0f, star_size, 0.0f);
        coot::Cartesian vec_z(0.0f, 0.0f, star_size);

        int col = 0;
        if (atom_colour_type == 0 && atom_colour_map_p && cra.chain) {
            col = atom_colour_map_p->index_for_chain(cra.chain->name);
        }
        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

        add_bond_half_gemmi(col, atom_pos + vec_x, atom_pos - vec_x, cc, model_number, -1, -1, true, true);
        add_bond_half_gemmi(col, atom_pos + vec_y, atom_pos - vec_y, cc, model_number, -1, -1, true, true);
        add_bond_half_gemmi(col, atom_pos + vec_z, atom_pos - vec_z, cc, model_number, -1, -1, true, true);

        bonds_added = 3;
    }

    return bonds_added;
}
