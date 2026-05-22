#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include "coot-utils/cablam-markup.hh"

namespace coot {

// Ported constructor for cablam_markup_t using gemmi types
// The original MMDB constructor was:
//   cablam_markup_t() { score = -1; residue = 0;}
// 
// We keep the same semantics but accept a gemmi::Residue* instead of mmdb::Residue*
class cablam_markup_t_gemmi {
public:
    clipper::Coord_orth O_prev_pos;
    clipper::Coord_orth O_this_pos;
    clipper::Coord_orth O_next_pos;
    clipper::Coord_orth CA_proj_point_prev;
    clipper::Coord_orth CA_proj_point_this;
    clipper::Coord_orth CA_proj_point_next;
    gemmi::Residue* residue;  // pointer to the central residue (gemmi version)
    double score;

    // Default constructor - matches original semantics
    cablam_markup_t_gemmi() { score = -1; residue = nullptr; }
    
    // Constructor that accepts gemmi residue pointer
    explicit cablam_markup_t_gemmi(gemmi::Residue* res) { 
        score = -1; 
        residue = res; 
    }
};

} // namespace coot