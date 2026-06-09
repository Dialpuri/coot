#pragma once
#include <string>
#include "geometry/protein-geometry.hh"

namespace coot {

inline std::string get_group_gemmi(const protein_geometry& geom, const std::string& res_name_in) {
    std::string res_name = res_name_in;

    // First loop: match by three_letter_code
    unsigned int s = geom.size();
    for (unsigned int i = 0; i < s; i++) {
        if (geom.three_letter_code(i) == res_name) {
            return geom.operator[](i).second.residue_info.group;
        }
    }

    // Second loop: match by comp_id (using public operator[] instead of private dict_res_restraints)
    for (unsigned int i = 0; i < s; i++) {
        if (geom.operator[](i).second.residue_info.comp_id == res_name) {
            return geom.operator[](i).second.residue_info.group;
        }
    }

    std::string msg = "WARNING:: get_group(): No dictionary group found for residue type :";
    msg += res_name;
    msg += ":";
    throw std::runtime_error(msg);
}

} // namespace coot
