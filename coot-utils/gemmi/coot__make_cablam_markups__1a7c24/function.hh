#pragma once
#include <vector>
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

struct cablam_markup_t {
  clipper::Coord_orth O_prev_pos;
  clipper::Coord_orth O_this_pos;
  clipper::Coord_orth O_next_pos;
  clipper::Coord_orth CA_proj_point_prev;
  clipper::Coord_orth CA_proj_point_this;
  clipper::Coord_orth CA_proj_point_next;
  gemmi::CRA cra;
  double score;
  
  cablam_markup_t() : score(-1) {}
};

std::vector<cablam_markup_t>
make_cablam_markups_gemmi(const std::vector<std::pair<residue_spec_t, double> > &residues,
                          gemmi::Structure &st);

} // namespace coot