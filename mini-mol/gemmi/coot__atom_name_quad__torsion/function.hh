#pragma once

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <clipper/clipper.h>

namespace coot {

struct atom_name_quad {
    std::string atom_name_[4];

    atom_name_quad(const std::string& a0, const std::string& a1,
                   const std::string& a2, const std::string& a3) {
        atom_name_[0] = a0;
        atom_name_[1] = a1;
        atom_name_[2] = a2;
        atom_name_[3] = a3;
    }

    double torsion_gemmi(const gemmi::Residue& residue) const {
        double r = -999.9;

        const gemmi::Atom* at_0 = nullptr;
        const gemmi::Atom* at_1 = nullptr;
        const gemmi::Atom* at_2 = nullptr;
        const gemmi::Atom* at_3 = nullptr;

        for (const auto& atom : residue.atoms) {
            if (atom.name == atom_name_[0]) at_0 = &atom;
            else if (atom.name == atom_name_[1]) at_1 = &atom;
            else if (atom.name == atom_name_[2]) at_2 = &atom;
            else if (atom.name == atom_name_[3]) at_3 = &atom;
        }

        if (at_0 && at_1 && at_2 && at_3) {
            clipper::Coord_orth pt_0(at_0->pos.x, at_0->pos.y, at_0->pos.z);
            clipper::Coord_orth pt_1(at_1->pos.x, at_1->pos.y, at_1->pos.z);
            clipper::Coord_orth pt_2(at_2->pos.x, at_2->pos.y, at_2->pos.z);
            clipper::Coord_orth pt_3(at_3->pos.x, at_3->pos.y, at_3->pos.z);
            double angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_0, pt_1, pt_2, pt_3));
            return angle;
        }
        return r;
    }
};

} // namespace coot