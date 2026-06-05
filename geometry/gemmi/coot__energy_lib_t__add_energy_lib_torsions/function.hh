#pragma once

#include <coot/geometry/energy-lib.hh>
#include <gemmi/cif.hpp>
#include <string>
#include <vector>

namespace coot {

inline void add_energy_lib_torsions_gemmi(energy_lib_t& energy_lib, gemmi::cif::Loop loop) {
    size_t num_rows = loop.length();

    for (size_t j = 0; j < num_rows; j++) {
        std::string atom_type_1;
        std::string atom_type_2;
        std::string atom_type_3;
        std::string atom_type_4;
        double constant = 0;
        double angle = 0;
        int period = 0;
        bool ok = true;

        // Get column indices
        int col_atom_type_1 = loop.find_tag("atom_type_1");
        int col_atom_type_2 = loop.find_tag("atom_type_2");
        int col_atom_type_3 = loop.find_tag("atom_type_3");
        int col_atom_type_4 = loop.find_tag("atom_type_4");
        int col_const       = loop.find_tag("const");
        int col_angle       = loop.find_tag("angle");
        int col_period      = loop.find_tag("period");

        // Check that all columns exist
        if (col_atom_type_1 < 0 || col_atom_type_2 < 0 || col_atom_type_3 < 0 ||
            col_atom_type_4 < 0 || col_const < 0 || col_angle < 0 || col_period < 0) {
            continue;
        }

        // Read string columns
        atom_type_1 = loop.val(j, static_cast<size_t>(col_atom_type_1));
        atom_type_2 = loop.val(j, static_cast<size_t>(col_atom_type_2));
        atom_type_3 = loop.val(j, static_cast<size_t>(col_atom_type_3));
        atom_type_4 = loop.val(j, static_cast<size_t>(col_atom_type_4));

        // Read const (real)
        try {
            constant = std::stod(loop.val(j, static_cast<size_t>(col_const)));
        } catch (...) {
            ok = false;
        }

        // Read angle (real)
        try {
            angle = std::stod(loop.val(j, static_cast<size_t>(col_angle)));
        } catch (...) {
            ok = false;
        }

        // Read period (integer)
        try {
            period = std::stoi(loop.val(j, static_cast<size_t>(col_period)));
        } catch (...) {
            ok = false;
        }

        if (ok) {
            energy_lib_torsion tors(atom_type_1, atom_type_2, atom_type_3, atom_type_4,
                                    constant, angle, period);
            energy_lib.torsions.push_back(tors);
        }
    }
}

} // namespace coot