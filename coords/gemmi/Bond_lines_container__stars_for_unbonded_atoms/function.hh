#pragma once

#include "coords/Bond_lines.hh"
#include <gemmi/model.hpp>
#include <unordered_set>

// Helper class that inherits from Bond_lines_container to expose protected members
class Bond_lines_container_gemmi : public Bond_lines_container {
public:
    using Bond_lines_container::bonds;

    void add_bond_gemmi(int colour, const coot::Cartesian& first, const coot::Cartesian& second,
                        graphics_line_t::cylinder_class_t cc,
                        int chain_a, int chain_b, int chain_c,
                        bool draw_a, bool draw_b,
                        bool draw_a_b, bool draw_b_a)
    {
        addBond(colour, first, second, cc, chain_a, chain_b, chain_c, draw_a, draw_b);
    }

    // Count total bonds across all bond groups
    int total_bond_count() const {
        int total = 0;
        for (const auto& bl : bonds) total += bl.size();
        return total;
    }
};

// gemmi port of Bond_lines_container::stars_for_unbonded_atoms(mmdb::Manager*, int)
// Instead of MMDB UDD handles, passes a set of pointers to atoms that should receive stars.
inline void stars_for_unbonded_atoms_gemmi(
    Bond_lines_container_gemmi& bonds,
    const gemmi::Model& model,
    const std::unordered_set<const gemmi::Atom*>& unbonded_atoms)
{
    float star_size = 0.2f;
    graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

    coot::Cartesian small_vec_x(star_size, 0.1f, 0.1f);
    coot::Cartesian small_vec_y(0.1f, star_size, 0.1f);
    coot::Cartesian small_vec_z(0.1f, 0.1f, star_size);

    int col = 2;

    for (const auto& chain : model.chains) {
        for (const auto& res : chain.residues) {
            for (const auto& atom : res.atoms) {
                if (unbonded_atoms.count(&atom)) {
                    coot::Cartesian atom_pos(atom.pos.x, atom.pos.y, atom.pos.z);
                    bonds.add_bond_gemmi(col, atom_pos + small_vec_x, atom_pos - small_vec_x, cc, -1, -1, -1, true, true, true, true);
                    bonds.add_bond_gemmi(col, atom_pos + small_vec_y, atom_pos - small_vec_y, cc, -1, -1, -1, true, true, true, true);
                    bonds.add_bond_gemmi(col, atom_pos + small_vec_z, atom_pos - small_vec_z, cc, -1, -1, -1, true, true, true, true);
                }
            }
        }
    }
}