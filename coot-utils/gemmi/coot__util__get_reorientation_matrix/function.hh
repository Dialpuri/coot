#pragma once

#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include <utility>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {
namespace util {

std::pair<bool, clipper::RTop_orth>
get_reorientation_matrix_gemmi(const gemmi::Residue* residue_current,
                               const gemmi::Residue* residue_next) {

  bool status = false;
  clipper::RTop_orth rtop;

  if (residue_current) {
    if (residue_next) {
      const gemmi::Atom* atoms_this[3] = {nullptr, nullptr, nullptr};
      const gemmi::Atom* atoms_next[3] = {nullptr, nullptr, nullptr};

      // Find N, CA, C in current residue (gemmi stores names unpadded)
      for (const auto& at : residue_current->atoms) {
        if (at.name == "N")  atoms_this[0] = &at;
        if (at.name == "CA") atoms_this[1] = &at;
        if (at.name == "C")  atoms_this[2] = &at;
      }

      if (atoms_this[0] && atoms_this[1] && atoms_this[2]) {
        // Find N, CA, C in next residue
        for (const auto& at : residue_next->atoms) {
          if (at.name == "N")  atoms_next[0] = &at;
          if (at.name == "CA") atoms_next[1] = &at;
          if (at.name == "C")  atoms_next[2] = &at;
        }

        if (atoms_next[0] && atoms_next[1] && atoms_next[2]) {
          std::vector<clipper::Coord_orth> this_pos;
          std::vector<clipper::Coord_orth> next_pos;
          this_pos.push_back(coot::co_gemmi(atoms_this[0]));
          this_pos.push_back(coot::co_gemmi(atoms_this[1]));
          this_pos.push_back(coot::co_gemmi(atoms_this[2]));
          next_pos.push_back(coot::co_gemmi(atoms_next[0]));
          next_pos.push_back(coot::co_gemmi(atoms_next[1]));
          next_pos.push_back(coot::co_gemmi(atoms_next[2]));
          clipper::RTop_orth lrtop(this_pos, next_pos);
          rtop = lrtop;
          status = true;
        }
      }
    }
  }

  return std::pair<bool, clipper::RTop_orth>(status, rtop);
}

} // namespace util
} // namespace coot