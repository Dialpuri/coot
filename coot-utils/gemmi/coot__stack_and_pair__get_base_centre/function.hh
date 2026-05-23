#pragma once

#include <utility>
#include <string>
#include <set>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>

namespace coot {
namespace stack_and_pair {

// Get the centre of base atoms for a nucleic acid residue
// Returns (success, centre) where centre is in clipper::Coord_orth coordinates
inline clipper::Coord_orth co_gemmi(const gemmi::Atom* at) {
    return clipper::Coord_orth(at->pos.x, at->pos.y, at->pos.z);
}

// Initialize base atom name set
inline std::set<std::string> init_base_atom_names() {
    std::set<std::string> names;
    // Exact copy from coot::stack_and_pair::init()
    // PDBv3 FIXME - note the spaces in atom names
    names.insert(" N1 "); names.insert(" C2 "); names.insert(" N3 ");
    names.insert(" C4 "); names.insert(" C5 "); names.insert(" C6 ");
    names.insert(" N7 "); names.insert(" C8 "); names.insert(" N9 ");
    return names;
}

// Get the centre of base atoms for a nucleic acid residue
// Returns (success, centre) where centre is in clipper::Coord_orth coordinates
inline std::pair<bool, clipper::Coord_orth>
get_base_centre_gemmi(const gemmi::Residue* residue) {
    std::pair<bool, clipper::Coord_orth> p(false, clipper::Coord_orth(0,0,0));
    
    static const std::set<std::string> base_atom_name_set = init_base_atom_names();
    
    unsigned int n_centres = 0;
    clipper::Coord_orth centre_sum(0.0, 0.0, 0.0);
    
    for (const gemmi::Atom& atom : residue->atoms) {
        if (base_atom_name_set.find(atom.name) != base_atom_name_set.end()) {
            centre_sum += co_gemmi(&atom);
            n_centres++;
        }
    }
    
    if (n_centres > 3) {
        p.first = true;
        double n_d = static_cast<double>(n_centres);
        double r_n_d = 1.0 / n_d;
        p.second = clipper::Coord_orth(centre_sum.x() * r_n_d,
                                       centre_sum.y() * r_n_d,
                                       centre_sum.z() * r_n_d);
    }
    return p;
}

} // namespace stack_and_pair
} // namespace coot