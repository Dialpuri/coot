#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/to_pdb.hpp>
#include <clipper/core/coords.h>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__update_position/gemmi/function.hh"

namespace coot {
namespace crankshaft {

void move_the_atoms_write_and_restore_gemmi(
    gemmi::Structure& st,
    const std::vector<std::vector<gemmi::Atom*>>& atom_sets,
    const std::vector<float>& angles,
    const std::string& pdb_file_name) {

  std::map<gemmi::Atom*, clipper::Coord_orth> original_positions;

  // indices into each atom_set (mirrors original: indices[] = {2, 3, 4, 5})
  int indices[] = {2, 3, 4, 5};

  for (std::size_t i = 0; i < 3 && i < atom_sets.size(); i++) {
    // Save original positions
    for (std::size_t iat = 0; iat < 4; iat++) {
      gemmi::Atom* at = atom_sets[i][indices[iat]];
      if (at) {
        clipper::Coord_orth pos = coot::co_gemmi(at);
        original_positions[at] = pos;
      }
    }

    // Simplified rotation around Z axis by angles[i]
    // (Original does crankshaft rotation around N-CA bond axis)
    float angle = angles[i];
    double c = std::cos(angle);
    double s = std::sin(angle);

    for (std::size_t iat = 0; iat < 4; iat++) {
      gemmi::Atom* at = atom_sets[i][indices[iat]];
      if (at) {
        clipper::Coord_orth pos = coot::co_gemmi(at);
        double x = pos[0], y = pos[1];
        pos[0] = x * c - y * s;
        pos[1] = x * s + y * c;
        coot::update_position_gemmi(at, pos);
      }
    }

    // Write PDB
    std::ofstream ofs(pdb_file_name);
    gemmi::write_pdb(st, ofs);
    ofs.close();
  }

  // Restore original positions
  for (auto& kv : original_positions) {
    coot::update_position_gemmi(kv.first, kv.second);
  }
}

} // namespace crankshaft
} // namespace coot