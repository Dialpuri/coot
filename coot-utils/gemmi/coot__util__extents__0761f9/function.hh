#pragma once
#include <string>
#include <utility>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <clipper/core/coords.h>

namespace coot { namespace util {

// Ported from MMDB version: calculates bounding box extents for selected atoms
// Instead of SelectionHandle, this takes a gemmi::Structure and optional chain filter
std::pair<clipper::Coord_orth, clipper::Coord_orth>
extents_gemmi(const gemmi::Structure &st,
              const std::string &chain_id = "") {
  float most_x = -99999.f;
  float most_y = -99999.f;
  float most_z = -99999.f;
  float least_x = 99999.f;
  float least_y = 99999.f;
  float least_z = 99999.f;
  
  bool use_chain_filter = !chain_id.empty() && chain_id != "*";
  
  for (const gemmi::Model& model : st.models) {
    for (const gemmi::Chain& chain : model.chains) {
      if (use_chain_filter && chain.name != chain_id)
        continue;
      for (const gemmi::Residue& residue : chain.residues) {
        for (const gemmi::Atom& atom : residue.atoms) {
          float x = atom.pos.x;
          float y = atom.pos.y;
          float z = atom.pos.z;
          
          if (x > most_x) most_x = x;
          if (y > most_y) most_y = y;
          if (z > most_z) most_z = z;
          
          if (x < least_x) least_x = x;
          if (y < least_y) least_y = y;
          if (z < least_z) least_z = z;
        }
      }
    }
  }
  
  clipper::Coord_orth p1(most_x, most_y, most_z);
  clipper::Coord_orth p2(least_x, least_y, least_z);
  
  return std::pair<clipper::Coord_orth, clipper::Coord_orth>(p2, p1);
}

}} // namespace coot::util