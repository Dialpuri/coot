#pragma once

#include <vector>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>
#include <stdexcept>
#include <string>

namespace coot {

class atom_index_quad {
public:
    int index1, index2, index3, index4;

    atom_index_quad(int i1, int i2, int i3, int i4)
        : index1(i1), index2(i2), index3(i3), index4(i4) {}

    double torsion_gemmi(const std::vector<gemmi::Atom>& atoms, int n_selected_atoms) const {
        double angle = 0;
        for (int i = 0; i < n_selected_atoms; i++) {
            bool good_indexing = false;
            if ((index1 >= 0) && (index1 < n_selected_atoms)) {
                if ((index2 >= 0) && (index2 < n_selected_atoms)) {
                    if ((index3 >= 0) && (index3 < n_selected_atoms)) {
                        if ((index4 >= 0) && (index4 < n_selected_atoms)) {
                            good_indexing = true;
                        }
                    }
                }
            }
            if (!good_indexing) {
                std::string mess = "bad atom indexing in atom_index_quad::torsion()";
                throw std::runtime_error(mess);
            }
            clipper::Coord_orth pt_1(atoms[index1].pos.x, atoms[index1].pos.y, atoms[index1].pos.z);
            clipper::Coord_orth pt_2(atoms[index2].pos.x, atoms[index2].pos.y, atoms[index2].pos.z);
            clipper::Coord_orth pt_3(atoms[index3].pos.x, atoms[index3].pos.y, atoms[index3].pos.z);
            clipper::Coord_orth pt_4(atoms[index4].pos.x, atoms[index4].pos.y, atoms[index4].pos.z);

            angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4));
        }
        return angle;
    }
};

} // namespace coot