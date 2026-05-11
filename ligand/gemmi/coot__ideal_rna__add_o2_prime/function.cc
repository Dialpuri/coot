#include "function.hh"
#include <gemmi/pdb.hpp>
#include <gemmi/elem.hpp>

namespace coot {
namespace ideal_rna {

void add_o2_prime_gemmi(gemmi::Residue* res) {
    if (!res)
        return;

    // Reference points in the ideal geometry (C1', C2', C3')
    std::vector<clipper::Coord_orth> mov_pts;
    mov_pts.push_back(clipper::Coord_orth(6.853, -5.219, 1.725));  // C1'
    mov_pts.push_back(clipper::Coord_orth(7.473, -4.888, 3.087));  // C2'
    mov_pts.push_back(clipper::Coord_orth(6.678, -5.807, 4.033));  // C3'

    // The coordinate that gets moved by the rtop returned from LSQing.
    clipper::Coord_orth o2p(8.870, -5.158, 3.018);

    // Find C1', C2', C3' atoms in the residue
    gemmi::Atom* c1p = nullptr;
    gemmi::Atom* c2p = nullptr;
    gemmi::Atom* c3p = nullptr;

    for (auto& atom : res->atoms) {
        if (atom.name == " C1'" || atom.name == " C1*")
            c1p = &atom;
        if (atom.name == " C2'" || atom.name == " C2*")
            c2p = &atom;
        if (atom.name == " C3'" || atom.name == " C3*")
            c3p = &atom;
    }

    if (c1p && c2p && c3p) {
        std::vector<clipper::Coord_orth> ref_pts;
        ref_pts.push_back(clipper::Coord_orth(c1p->pos.x, c1p->pos.y, c1p->pos.z));
        ref_pts.push_back(clipper::Coord_orth(c2p->pos.x, c2p->pos.y, c2p->pos.z));
        ref_pts.push_back(clipper::Coord_orth(c3p->pos.x, c3p->pos.y, c3p->pos.z));

        clipper::RTop_orth rtop(mov_pts, ref_pts);
        clipper::Coord_orth pos = o2p.transform(rtop);

        gemmi::Atom new_atom;
        new_atom.name = " O2*";
        new_atom.element = gemmi::Element("O");
        new_atom.pos = gemmi::Position(pos.x(), pos.y(), pos.z());
        new_atom.occ = 1.0f;
        new_atom.b_iso = 30.0f;
        res->atoms.push_back(new_atom);
    }
}

}
}