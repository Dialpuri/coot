#pragma once
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include "coot-utils/coot-coord-utils.hh"

namespace coot {

// Port of position_residue_by_internal_coordinates::get_atom to gemmi
// Returns nullptr if atom not found
gemmi::Atom *
position_residue_by_internal_coordinates_get_atom_gemmi(gemmi::Residue *res_1,
                                                        gemmi::Residue *res_2,
                                                        const atom_name_quad &quad,
                                                        int atom_index) {

    gemmi::Atom *at = nullptr;

    std::string atom_name;

    if (atom_name.length() > 0) {

        gemmi::Residue *residue_p = res_1;
        std::string torsion_atom_name = quad.atom_name(atom_index);
        if (quad.atom_residue_index[atom_index] == 2)
            residue_p = res_2;
        
        for (gemmi::Atom& atom : residue_p->atoms) {
            if (atom.name == torsion_atom_name) {
                at = &atom;
                break;
            }
        }
    }
    return at;
}

} // namespace coot