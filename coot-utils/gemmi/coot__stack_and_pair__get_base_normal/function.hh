#pragma once

#include <string>
#include <vector>
#include <set>
#include <utility>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "coot-utils/coot-coord-utils.hh"

namespace coot {
namespace stack_and_pair {

static const std::set<std::string> base_atom_name_set = {
    " N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 ", " N7 ", " C8 ", " N9 ",
    " O2 ", " O4 ", " O6 ",
    " C4'", " C1'"
};

std::pair<bool, clipper::Coord_orth>
get_base_normal_gemmi(const gemmi::Residue* residue_p) {

  std::pair<bool, clipper::Coord_orth> r(false, clipper::Coord_orth(0, 0, 0));

  std::vector<clipper::Coord_orth> v;
  for (const gemmi::Atom& at : residue_p->atoms) {
    std::string atom_name(at.name);
    if (base_atom_name_set.find(atom_name) != base_atom_name_set.end()) {
      v.push_back(coot::co_gemmi(&at));
    }
  }

  if (v.size() > 2) {
    coot::lsq_plane_info_t lsq(v);
    r.first = true;
    r.second = clipper::Coord_orth(lsq.normal().unit());
  }
  return r;
}

} // namespace stack_and_pair
} // namespace coot