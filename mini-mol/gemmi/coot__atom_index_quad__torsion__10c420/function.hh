#pragma once

#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <string>

namespace coot {

struct atom_index_quad {
    int index1, index2, index3, index4;

    atom_index_quad(int i1, int i2, int i3, int i4)
        : index1(i1), index2(i2), index3(i3), index4(i4) {}

    double torsion_gemmi(const gemmi::Residue& residue) const {
        if (index1 < 0 || index1 >= static_cast<int>(residue.atoms.size()) ||
            index2 < 0 || index2 >= static_cast<int>(residue.atoms.size()) ||
            index3 < 0 || index3 >= static_cast<int>(residue.atoms.size()) ||
            index4 < 0 || index4 >= static_cast<int>(residue.atoms.size())) {
            throw std::runtime_error("bad atom indexing in atom_index_quad::torsion()");
        }

        const auto& pt_1 = residue.atoms[index1].pos;
        const auto& pt_2 = residue.atoms[index2].pos;
        const auto& pt_3 = residue.atoms[index3].pos;
        const auto& pt_4 = residue.atoms[index4].pos;

        return gemmi::calculate_dihedral(pt_1, pt_2, pt_3, pt_4) * 180.0 / M_PI;
    }
};

} // namespace coot