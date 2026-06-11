#pragma once

#include <string>
#include <vector>
#include <map>

#include "geometry/protein-geometry.hh"

// Pad an atom name to 4-char PDB column format (identical logic to coot::atom_id_mmdb_expand).
inline std::string atom_id_mmdb_expand(const std::string& atomname) {
    std::string r;
    int ilen = static_cast<int>(atomname.length());

    if (ilen == 4) return atomname;

    if (ilen == 1) {
        r = " ";
        r += atomname;
        r += "  ";
    } else {
        if (ilen == 2) {
            static const char* metals[] = {
                "MG","NA","LI","AL","SI","CL","SC","TI","CR","MN","FE","CO",
                "NI","CU","ZN","GA","AS","SE","BR","RB","SR","RE","OS","IR",
                "PT","AU","HG","PB","BI"
            };
            bool done = false;
            for (const char* m : metals) {
                if (atomname == m) { done = true; break; }
            }
            if (done) {
                r += atomname;
                r += "  ";
            } else {
                r = " ";
                r += atomname;
                r += " ";
            }
        } else {
            if (ilen == 3) {
                r = " ";
                r += atomname;
            } else {
                r = atomname;
            }
        }
    }
    return r;
}

namespace coot {

inline void add_chem_mod_angle_gemmi(
    const std::vector<std::string>& tags,
    const std::vector<std::vector<std::string>>& rows,
    std::map<std::string, coot::chem_mod>& mods) {

    auto find_col = [&tags](const std::string& tag) -> int {
        for (int i = 0; i < static_cast<int>(tags.size()); ++i) {
            if (tags[i] == tag) return i;
        }
        return -1;
    };

    int n_rows = static_cast<int>(rows.size());
    int col_mod_id        = find_col("mod_id");
    int col_function      = find_col("function");
    int col_atom_id_1     = find_col("atom_id_1");
    int col_atom_id_2     = find_col("atom_id_2");
    int col_atom_id_3     = find_col("atom_id_3");
    int col_new_value     = find_col("new_value_angle");
    int col_new_value_esd = find_col("new_value_angle_esd");

    for (int j = 0; j < n_rows; j++) {
        int ierr_tot = 0;

        std::string mod_id;
        std::string function;
        std::string atom_id_1;
        std::string atom_id_2;
        std::string atom_id_3;
        double new_value_angle;
        double new_value_angle_esd;

        // mod_id
        if (col_mod_id < 0 || static_cast<size_t>(col_mod_id) >= rows[j].size()) {
            ierr_tot++;
        } else {
            mod_id = rows[j][col_mod_id];
        }

        // function
        if (col_function < 0 || static_cast<size_t>(col_function) >= rows[j].size()) {
            ierr_tot++;
        } else {
            function = rows[j][col_function];
        }

        // atom_id_1
        if (col_atom_id_1 < 0 || static_cast<size_t>(col_atom_id_1) >= rows[j].size()) {
            ierr_tot++;
        } else {
            atom_id_1 = rows[j][col_atom_id_1];
        }

        // atom_id_2
        if (col_atom_id_2 < 0 || static_cast<size_t>(col_atom_id_2) >= rows[j].size()) {
            ierr_tot++;
        } else {
            atom_id_2 = rows[j][col_atom_id_2];
        }

        // atom_id_3
        if (col_atom_id_3 < 0 || static_cast<size_t>(col_atom_id_3) >= rows[j].size()) {
            ierr_tot++;
        } else {
            atom_id_3 = rows[j][col_atom_id_3];
        }

        // new_value_angle
        if (col_new_value < 0 || static_cast<size_t>(col_new_value) >= rows[j].size()) {
            ierr_tot++;
        } else {
            std::string sv = rows[j][col_new_value];
            if (sv.empty() || sv == ".") {
                ierr_tot++;
            } else {
                new_value_angle = std::stod(sv);
            }
        }

        // new_value_angle_esd
        if (col_new_value_esd < 0 || static_cast<size_t>(col_new_value_esd) >= rows[j].size()) {
            ierr_tot++;
        } else {
            std::string sv = rows[j][col_new_value_esd];
            if (sv.empty() || sv == ".") {
                ierr_tot++;
            } else {
                new_value_angle_esd = std::stod(sv);
            }
        }

        if (ierr_tot == 0) {
            coot::chem_mod_angle cma(function,
                                     atom_id_mmdb_expand(atom_id_1),
                                     atom_id_mmdb_expand(atom_id_2),
                                     atom_id_mmdb_expand(atom_id_3),
                                     new_value_angle,
                                     new_value_angle_esd);
            mods[mod_id].add_mod_angle(cma);
        }
    }
}

} // namespace coot
