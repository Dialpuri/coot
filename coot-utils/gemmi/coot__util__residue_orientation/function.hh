#pragma once
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/clipper.h>
#include <vector>
#include <string>
#include <iostream>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__is_main_chain_p/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__average_position/gemmi/function.hh"

namespace coot {
namespace util {

inline clipper::Mat33<double>
residue_orientation_gemmi(const gemmi::Residue& residue, const clipper::Mat33<double>& orientation_in) {
    std::vector<clipper::Coord_orth> pts;
    clipper::Mat33<double> r(0,0,0,0,0,0,0,0,0);  // Initialize to zero matrix
    clipper::Coord_orth n_vec(0,0,1);
    const gemmi::Atom* ca = nullptr;
    const gemmi::Atom* n = nullptr;
    
    for (const gemmi::Atom& atom : residue.atoms) {
        if (!coot::is_main_chain_p_gemmi(atom, residue)) {
            clipper::Coord_orth coord(atom.pos.x, atom.pos.y, atom.pos.z);
            pts.push_back(coord);
        } else {
            if (atom.name == "CA") {
                ca = &atom;
            }
            if (atom.name == "N") {
                n = &atom;
            }
        } 
    }
    
    if (pts.size() > 0) {
        if (ca) { 
            clipper::Coord_orth ca_pos(ca->pos.x, ca->pos.y, ca->pos.z);
            std::vector<clipper::Coord_orth> pts_copy = pts;
            clipper::Coord_orth average_pos = coot::util::average_position_gemmi(pts_copy);
            clipper::Coord_orth u((average_pos - ca_pos).unit());

            // reset n_vec to something sensible, if we have the CA and N.
            if (ca && n) {
                clipper::Coord_orth n_pos(n->pos.x, n->pos.y, n->pos.z);
                n_vec = n_pos - ca_pos;
            }
        
            // now make a mat...
            clipper::Coord_orth n_vec_unit(n_vec.unit());
            
            clipper::Coord_orth p1(clipper::Coord_orth::cross(n_vec_unit, u).unit());
            clipper::Coord_orth p2(clipper::Coord_orth::cross(p1, u).unit());
            clipper::Coord_orth p3 = u;

            r = clipper::Mat33<double>(p1.x(), p1.y(), p1.z(),
                                       p2.x(), p2.y(), p2.z(),
                                       p3.x(), p3.y(), p3.z());
        }
    }
    return r;
}

}
}