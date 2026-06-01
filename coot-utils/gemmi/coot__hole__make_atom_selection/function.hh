#pragma once

#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "clipper/core/coords.h"
#include <vector>
#include <string>
#include <set>

namespace coot {
namespace hole {

/*  Return CRAs (Chain*, Residue*, Atom*) whose element is one of
    H/C/O/N/S, whose residue name is not HOH, and whose position
    is within *radius* of *pt*.

    This is the gemmi-equivalent of the original MMDB function
    coot::hole::make_atom_selection(), which used an MMDB selection
    handle to collect the same set of atoms.
*/
inline std::vector<gemmi::CRA>
make_atom_selection_gemmi(const gemmi::Model& model,
                          const clipper::Coord_orth& pt,
                          double radius) {
  std::vector<gemmi::CRA> result;

  double radius_sq = radius * radius;
  double px = pt.x(), py = pt.y(), pz = pt.z();

  // Allowed elements as strings
  std::set<std::string> allowed_elems = {"H", "C", "O", "N", "S"};

  for (const gemmi::Chain& chain : model.chains) {
    for (const gemmi::Residue& res : chain.residues) {
      // Skip HOH residues
      if (res.name == "HOH") continue;
      for (const gemmi::Atom& atom : res.atoms) {
        std::string elem_name = atom.element.name();
        if (allowed_elems.find(elem_name) == allowed_elems.end()) continue;

        double dx = atom.pos.x - px;
        double dy = atom.pos.y - py;
        double dz = atom.pos.z - pz;
        if (dx*dx + dy*dy + dz*dz <= radius_sq) {
          result.push_back({const_cast<gemmi::Chain*>(&chain),
                            const_cast<gemmi::Residue*>(&res),
                            const_cast<gemmi::Atom*>(&atom)});
        }
      }
    }
  }

  return result;
}

} // namespace hole
} // namespace coot