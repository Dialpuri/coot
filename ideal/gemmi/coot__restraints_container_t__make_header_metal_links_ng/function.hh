#pragma once

#include <string>
#include <vector>
#include "geometry/protein-geometry.hh"

// Replacement for mmdb::Link — holds the data from a structural link
struct link_info_t {
  std::string chainID1;
  int         seqNum1;
  std::string insCode1;
  std::string atName1;
  std::string resName1;
  std::string aloc1;

  std::string chainID2;
  int         seqNum2;
  std::string insCode2;
  std::string atName2;
  std::string resName2;
  std::string aloc2;
};

namespace coot {

// Port of coot::restraints_container_t::make_header_metal_links_ng
// Instead of mutating a restraints_container_t, returns the number of
// metal-link bonds that would have been added.  Signature mirrors the
// original logic: iterate links, look up element types, check metal maps.
int make_header_metal_links_ng_gemmi(
    const std::vector<link_info_t> &links,
    const coot::protein_geometry &geom);

} // namespace coot