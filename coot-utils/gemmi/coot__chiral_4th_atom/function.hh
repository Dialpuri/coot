#pragma once
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

// Port of coot::chiral_4th_atom from MMDB to gemmi
// Finds the closest atom to the center atom in a chiral center
// Returns nullptr if no atom is found within the distance threshold
inline gemmi::Atom*
chiral_4th_atom_gemmi(gemmi::Residue* residue, const gemmi::Atom* at_centre,
                      const gemmi::Atom* at_1, const gemmi::Atom* at_2,
                      const gemmi::Atom* at_3) {
    gemmi::Atom* rat = nullptr;
    double d_crit = std::sqrt(1.7);
    double d_sqrd = d_crit * d_crit; // tracks "best/closest"
    
    clipper::Coord_orth p_c = co_gemmi(at_centre);
    
    for (const gemmi::Atom& at : residue->atoms) {
        if (&at != at_centre && &at != at_1 && &at != at_2 && &at != at_3) {
            clipper::Coord_orth pt = co_gemmi(&at);
            double d2 = (p_c - pt).lengthsq();
            if (d2 < d_sqrd) {
                rat = const_cast<gemmi::Atom*>(&at);
                d_sqrd = d2;
            }
        }
    }
    
    return rat;
}

} // namespace coot