#pragma once

#include "coot/mini-mol/mini-mol.hh"
#include <gemmi/model.hpp>
#include <string>

namespace coot { namespace minimol {

inline atom atom_gemmi(gemmi::Atom const& at) {
  atom a;
  // Reconstruct the 4-char PDB atom name from gemmi's trimmed name.
  // Convention: letter-first names get a leading space (columns 14-16 in PDB),
  // digit-first names fill columns 13-16 directly.
  {
    std::string nm = at.name;
    if (nm.empty()) {
      a.name = "    ";
    } else if (std::isalpha(static_cast<unsigned char>(nm[0]))) {
      std::string padded = " " + nm;
      while (padded.size() < 4) padded += ' ';
      a.name = padded;
    } else {
      // digit-first
      std::string padded = nm;
      while (padded.size() < 4) padded += ' ';
      a.name = padded;
    }
  }

  // gemmi element.name() returns unpadded (e.g. "C"); MMDB stored left-padded to 2 chars (e.g. " C")
  std::string elem = at.element.name();
  while (elem.size() < 2)
    elem.insert(elem.begin(), ' ');
  a.element = elem;
  a.pos = clipper::Coord_orth(at.pos.x, at.pos.y, at.pos.z);
  // gemmi altloc is char; '\0' or ' ' means no alt-loc, same as MMDB ""
  if (at.altloc == '\0' || at.altloc == ' ')
    a.altLoc = "";
  else
    a.altLoc = std::string(1, at.altloc);
  a.occupancy         = at.occ;
  a.temperature_factor = at.b_iso;
  a.int_user_data     = -1;
  return a;
}

}} // namespace coot::minimol
