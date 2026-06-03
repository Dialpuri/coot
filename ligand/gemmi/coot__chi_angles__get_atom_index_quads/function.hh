#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>

namespace coot {

struct atom_name_quad {
    std::string n[4];
    atom_name_quad(const std::string& n0, const std::string& n1,
                   const std::string& n2, const std::string& n3)
        : n{n0, n1, n2, n3} {}
    std::string atom_name(int i) const {
        if (i < 0 || i > 3)
            throw std::runtime_error("out of bounds index on atom_name_quad::atom_name()");
        return n[i];
    }
};

struct atom_index_quad {
    int index1;
    int index2;
    int index3;
    int index4;
    atom_index_quad() : index1(-1), index2(-1), index3(-1), index4(-1) {}
    atom_index_quad(int i1, int i2, int i3, int i4)
        : index1(i1), index2(i2), index3(i3), index4(i4) {}
};

namespace chi_angles {

inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<coot::atom_index_quad>
get_atom_index_quads_gemmi(
    const std::vector<coot::atom_name_quad> &atom_name_quads,
    const std::vector<std::string> &atom_names) {

    std::vector<coot::atom_index_quad> index_quads;
    for (size_t iquad = 0; iquad < atom_name_quads.size(); iquad++) {
        int index_1 = -1;
        int index_2 = -1;
        int index_3 = -1;
        int index_4 = -1;
        for (int iat = 0; iat < static_cast<int>(atom_names.size()); iat++) {
            std::string t_atomname = trim(atom_names[iat]);
            for (int j = 0; j < 4; j++) {
                std::string t_query = trim(atom_name_quads[iquad].atom_name(j));
                if (t_atomname == t_query) {
                    if (j == 0) index_1 = iat;
                    else if (j == 1) index_2 = iat;
                    else if (j == 2) index_3 = iat;
                    else index_4 = iat;
                }
            }
        }
        if ((index_1 != -1) && (index_2 != -1) && (index_3 != -1) && (index_4 != -1)) {
            coot::atom_index_quad iq(index_1, index_2, index_3, index_4);
            index_quads.push_back(iq);
        }
    }
    return index_quads;
}

} // namespace chi_angles
} // namespace coot