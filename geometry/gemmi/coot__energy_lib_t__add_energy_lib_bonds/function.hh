#pragma once

#include "coot/geometry/energy-lib.hh"
#include <gemmi/cif.hpp>
#include <string>
#include <cstdlib>

namespace coot {

void add_energy_lib_bonds_gemmi(energy_lib_t& energy_lib, const gemmi::cif::Loop& loop) {

    int col_atom1 = loop.find_tag("atom_type_1");
    int col_atom2 = loop.find_tag("atom_type_2");
    int col_type  = loop.find_tag("type");
    int col_const = loop.find_tag("const");
    int col_len   = loop.find_tag("length");
    int col_esd   = loop.find_tag("value_esd");

    if (col_atom1 < 0 || col_atom2 < 0 || col_type < 0) return;

    size_t n = loop.length();
    for (size_t j = 0; j < n; ++j) {

        std::string atom_type_1 = loop.val(j, col_atom1);
        std::string atom_type_2 = loop.val(j, col_atom2);
        std::string type        = loop.val(j, col_type);

        double spring_const = 420.0;
        double length       = 0.0;
        double value_esd    = 0.02;
        int    ierr_tot     = 0;

        // Check string fields for emptiness (matches original GetString error behavior)
        if (atom_type_1.empty()) ++ierr_tot;
        if (atom_type_2.empty()) ++ierr_tot;
        if (type.empty())        ++ierr_tot;

        // GetReal("const"): on error, default to 420; does NOT add to ierr_tot
        if (col_const >= 0) {
            const std::string& c_text = loop.val(j, col_const);
            if (!c_text.empty()) {
                char* end = nullptr;
                double val = std::strtod(c_text.c_str(), &end);
                if (end != c_text.c_str()) {
                    spring_const = val;
                } else {
                    spring_const = 420.0;
                }
            } else {
                spring_const = 420.0;
            }
        }

        // GetReal("length"): error DOES add to ierr_tot
        if (col_len >= 0) {
            const std::string& l_text = loop.val(j, col_len);
            if (!l_text.empty()) {
                char* end = nullptr;
                double val = std::strtod(l_text.c_str(), &end);
                if (end != l_text.c_str()) {
                    length = val;
                } else {
                    ++ierr_tot;
                }
            } else {
                ++ierr_tot;
            }
        } else {
            ++ierr_tot;
        }

        // GetReal("value_esd"): on error, default to 0.02; does NOT add to ierr_tot
        if (col_esd >= 0) {
            const std::string& e_text = loop.val(j, col_esd);
            if (!e_text.empty()) {
                char* end = nullptr;
                double val = std::strtod(e_text.c_str(), &end);
                if (end != e_text.c_str()) {
                    value_esd = val;
                } else {
                    value_esd = 0.02;
                }
            } else {
                value_esd = 0.02;
            }
        }

        if (ierr_tot == 0) {
            // one bond (ST-OS) is delocal, not deloc
            if (type == "delocal") type = "deloc";
            coot::energy_lib_bond bond(atom_type_1, atom_type_2, type, spring_const,
                                       length, value_esd);
            energy_lib.add_energy_lib_bond(bond);
        }
    }
}

} // namespace coot