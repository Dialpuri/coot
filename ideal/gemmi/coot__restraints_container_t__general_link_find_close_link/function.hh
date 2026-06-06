#pragma once

#include <gemmi/model.hpp>
#include <utility>
#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__restraints_container_t__general_link_find_close_link_inner/gemmi/function.hh"

namespace coot {

// gemmi port of restraints_container_t::general_link_find_close_link
// Calls general_link_find_close_link_inner_gemmi twice (forward and reversed)
// and returns the first non-empty link name found, with the order_switch flag.
inline std::pair<std::string, bool>
general_link_find_close_link_gemmi(
    const std::vector<chem_link> &li,
    const gemmi::Residue *r1,
    const gemmi::Residue *r2,
    bool order_switch_flag,
    const protein_geometry &geom) {

  std::pair<std::string, bool> r("", order_switch_flag);
  std::string rs = general_link_find_close_link_inner_gemmi(li, r1, r2, order_switch_flag, geom);
  if (rs != "") {
    r.first = rs;
  } else {
    rs = general_link_find_close_link_inner_gemmi(li, r2, r1, order_switch_flag, geom);
    if (rs != "") {
      r.first = rs;
      r.second = true;
    }
  }
  return r;
}

} // namespace coot